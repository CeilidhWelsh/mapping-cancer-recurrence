# import packages and libraries
import os
import skrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from skrt import Patient, ROI, StructureSet, Image
from skrt.registration import get_default_pfiles, Registration, set_elastix_dir
from skrt.better_viewer import BetterViewer
from import_analysis import ImportPatient
from skrt.simulation import SyntheticImage


# additional functions for saving structure set 
def save_structure_set(structure_set, path, rewrite=True):
    """ check if the directory to save the structure set exists, if it does and rewrite is True, then write structure set, otherwise create directory and write the structure set """
    # check if the directory exists already
    # if the directory exists, skip
    # if not make the directory and write it

    if os.path.exists(path):
        if rewrite == True:
            structure_set.write(outdir=path)
    else:
        os.mkdir(path)
        structure_set.write(outdir=path)


# Patient Classes
class Patient_Info():
    """
    A class representing patient information and related operations.

    Parameters:
    - p: Patient object containing relevant patient CT, RTSTRUCT and RTDOSE information.
    - trial: String representing the trial name (i.e. 'import', 'voxtox').
    - recurrence_dict: Dictionary where the key is the new name for the recurrence, and the value is a list of possible names for
                        recurrence structure as labelled by the clinician.
    - ctv_names: Dictionary where key is the new name for CTV, and value is a list of possible names for Clinical Target Volumes.
    - comparison_struct_dict: Dictionary where keys are new names for structures used to compare registration performance, and values
                                are lists of possible names as labelled by the clinican

    Attributes:
    - p: Patient object containing relevant information.
    - trial: String representing the trial type ('import' or other).
    - fixed_image: String representing where the fixed image type is 'planning' or 'fixed'
    - recurrence_dict: Dictionary containing naming scheme for recurrence volume.
    - ctv_names: Dictionary  containing naming scheme for Clinical Target Volumes (CTVs).
    - comparison_struct_dict: Dictionary containing naming scheme for comparison structures.
    - fixed, moving: Image objects obtained from patient studies.
    - adjusted_fixed, adjusted_moving: Image objects obtained from patient studies with adjusted voxel sizing
    - planning_ss, relapse_ss, ctv_ss, comparison_ss: StructureSet objects obtained from patient studies.

    Methods:
    - get_images(): Obtains fixed and moving images for the patient, adjusting voxel sizes if needed.
    - get_structure_sets(): Obtains structure sets for planning, relapse, CTV, and comparison based on trial.
    """

    def __init__(self, p, trial, recurrence_dict, ctv_names, comparison_struct_dict):
        """ Initializes Patient_Info object with relevant information.

        Parameters:
        - p: Patient object containing relevant information.
        - trial: String representing the trial type ('import' or other).
        - recurrence_dict: Dictionary containing naming scheme for recurrence volume.
        - ctv_names: Dictionary  containing naming scheme for Clinical Target Volumes (CTVs).
        - comparison_struct_dict: Dictionary containing naming scheme for comparison structures.
        """

        # initalise variables
        self.p = p
        self.trial = trial
        self.fixed_image = 'planning'
        self.recurrence_dict = recurrence_dict
        self.ctv_names = ctv_names
        self.comparison_struct_dict = comparison_struct_dict

        # run methods
        self.fixed, self.moving, self.adjusted_fixed, self.adjusted_moving = self.get_images()
        self.planning_ss, self.relapse_ss, self.ctv_ss, self.comparison_ss = self.get_structure_sets()
        self.relapse = self.get_relapse_struct()
        self.dose = self.get_dose()

    # methods for obtaining patient images and patient structure sets

    def get_images(self):
        """Obtains fixed and moving images for the patient, checking for CT or MR format and adjusting voxel sizes if needed."""

        # logic for obtaining the images based on trial name
        if 'import' in self.trial:
            self.fixed = self.p.get_ct_plan()
            self.moving = self.p.get_ct_relapse()

        else:
            self.fixed = self.p.get_studies()[0].ct_images[-1]
            mr_structure_sets = self.p.studies[-1].get_structure_sets("mr")

            if mr_structure_sets:
                self.moving = self.p.get_studies()[-1].mr_images[-1]
            else:
                self.moving = self.p.get_studies()[-1].ct_images[-1]

        # match voxel size of images between the fixed and moving image
        self.adjusted_fixed, self.adjusted_moving = skrt.image.match_image_voxel_sizes(
            self.fixed, self.moving, voxel_size=(2, 2, 2))

        return self.fixed, self.moving, self.adjusted_fixed, self.adjusted_moving

    def get_structure_sets(self):
        """Obtains structure sets for planning, relapse, CTV, and comparison dependent on trial."""

        # logic for obtaining the planning structure set, comparison structure set and relapse structure set
        if 'import' in self.trial:
            self.planning_ss = self.p.ss_clinical
            self.relapse_ss = self.p.ss_relapse
            recurrence_ss = self.p.ss_recurrence.filtered_copy(
                self.recurrence_dict, keep_renamed_only=True)
            self.relapse_ss.add_roi(recurrence_ss.get_roi('tumour_volume'))
            self.comparison_ss = self.p.ss_plan
            ctv_ss = self.p.ss_clinical

        else:
            self.planning_structure_sets = self.p.studies[0].ct_structure_sets

            mr_structure_sets = self.p.studies[-1].get_structure_sets("mr")
            if mr_structure_sets:
                self.relapse_structure_sets = self.p.studies[-1].mr_structure_sets[-1]
            else:
                self.relapse_structure_sets = self.p.studies[-1].ct_structure_sets[-1]

            for ss in self.planning_structure_sets:
                if len(ss.filtered_copy(self.comparison_struct_dict, keep_renamed_only=True).get_roi_names()) != 0:
                    self.comparison_ss = ss

                if len(ss.filtered_copy({'CTV': ['*CTV*', '*ctv*', '*CTV_*', '*cvt*']}, keep_renamed_only=True).get_roi_names()) != 0:
                    self.planning_ss = ss
                    ctv_ss = ss.filtered_copy(
                        to_keep=['*CTV*', '*ctv*', '*CTV_*', '*cvt*'])

                if len(ss.filtered_copy({'left parotid': ['left parotid']}, keep_renamed_only=True).get_roi_names()) != 0:
                    self.clinical_ss = ss
                else:
                    self.clinical_ss = {}

            # obtain relapse structure set
            self.relapse_ss = self.relapse_structure_sets.filtered_copy(
                self.recurrence_dict, keep_renamed_only=True)

        # relabel structure sets
        self.planning_ss = self.planning_ss.filtered_copy(self.ctv_names)
        self.comparison_ss = self.comparison_ss.filtered_copy(
            self.comparison_struct_dict, keep_renamed_only=True)
        self.ctv_ss = ctv_ss.filtered_copy(self.ctv_names)

        # define associated images for each structure set
        self.ctv_ss.set_image(self.adjusted_fixed)
        self.comparison_ss.set_image(self.adjusted_fixed)
        self.planning_ss.set_image(self.adjusted_fixed)
        self.relapse_ss.set_image(self.adjusted_moving)

        return self.planning_ss, self.relapse_ss, self.ctv_ss, self.comparison_ss

    def get_relapse_struct(self):
        """Obtains the relapse structure from the relapse structure set for the trial.

        Returns:
        --------
        - relapse: ROI object representing the relapse structure relabelled 'tumour_volume'.
        """
        if 'import' in self.trial:
            self.relapse = self.p.ss_recurrence
        else:
            self.relapse = self.relapse_ss.get_roi('tumour_volume')

        return self.relapse

    def get_dose(self):
        """Obtains the dose information based on the trial.

        Returns:
        --------
        - dose: Dose object representing the patient's RTDOSE information.
        """

        if 'import' in self.trial:
            self.dose = self.p.get_dose_sum()
        elif self.trial == 'voxtox':
            self.dose = self.p.studies[0].la3_doses[0]
        else:
            self.dose = self.p.studies[0].ct_doses[0]

        self.dose.match_size(self.adjusted_fixed)

        return self.dose


# Registration Class
class Registration_Step(Patient_Info):
    """ A class representing a registration step for a patient, inheriting from the Patient_Info class.

    Parameters:
    - p: Patient object containing relevant information.
    - trial: String representing the trial type ('import' or other).
    - recurrence_dict:  Dictionary containing naming scheme for recurrence volume.
    - ctv_names:  Dictionary containing naming scheme for Clinical Target Volumes (CTVs).
    - comparison_struct_dict:  Dictionary containing naming scheme for comparison structure information.
    - pfile_dir: Path to the directory containing registration parameter files.
    - results_dir: Path to the directory for storing registration results.

    Attributes:
    - pfile_dir: Path to the directory containing registration parameter files.
    - results_dir: Path to the directory for storing registration results.

    Methods:
    - perform_reg(): Performs registration for a patient inheriting from the Patient_Info class.
    """

    def __init__(self, p, trial, recurrence_dict, ctv_names, comparison_struct_dict, pfile_dir, results_dir):
        super().__init__(p, trial, recurrence_dict, ctv_names, comparison_struct_dict)
        self.pfile_dir = pfile_dir
        self.results_dir = results_dir

    def perform_reg(self):
        """Performs registration for a patient inherting from the Patient_Info class"""
        #print('check', list(self.comparison_struct_dict.keys())[0])
        # instantiate reg as Registration Object
        print('checking in function:', self.results_dir, self.pfile_dir)

        reg = Registration(
            Path(f"{self.results_dir}"),
            fixed=self.adjusted_fixed,
            moving=self.adjusted_moving,
            # fixed_mask=p.get_mask_plan(),
            initial_alignment=list(self.comparison_struct_dict.keys())[
                0],  # this could use fixing
            pfiles={
                "translation": f'{self.pfile_dir}/MI_Translation.txt',
                "rigid": f'{self.pfile_dir}/MI_Rigid.txt',
                "bspline": f'{self.pfile_dir}/MI_BSpline30.txt'
            },
            overwrite=False,
            capture_output=True,
            # log_level = 'DEBUG'
        )

        # perform registration
        reg.register()

        # visualise registration results for b-spline step and save image to registration results directory
        for step in reg.steps:
            if step == 'bspline':
                print('> Registration Step Complete')
                reg.view_result(
                    step=step, save_as=f'{self.results_dir}/{step}/reg_image_{step}.png', zoom=True)


# Transform Structures
class Transform_Structures(Patient_Info):
    def __init__(self, p, trial, recurrence_dict, ctv_names, comparison_struct_dict, results_dir):
        super().__init__(p, trial, recurrence_dict, ctv_names, comparison_struct_dict)
        self.p = p
        self.results_dir = results_dir
        self.structure_path = f'{self.results_dir}/structures'
        self.ctv_names = ctv_names

    def write_structures(self):
        """ write structure sets for each patient to the correpsonding directory.

        Returns
        -------
        ROI object
            - ctv: an ROI object for the patient ctv/tumour-bed structure obtained from the planning structure set
        """

        # obtain ctv structure
        ctv = self.ctv_ss.get_roi(list(self.ctv_names.keys())[0])
        ctv.set_image(self.adjusted_fixed)

        # check structure path directory exists and write planning structures and ctv
        if not os.path.exists(f'{self.structure_path}/planning'):
            os.makedirs(f'{self.structure_path}/planning')
            self.planning_ss.write(outdir=f'{self.structure_path}/planning')
            ctv.write(outname='tumour_volume.nii',
                      outdir=f'{self.structure_path}/planning')

        # check structure path directory exists and write all relapse structures
        if not os.path.exists(f'{self.structure_path}/relapse'):
            os.makedirs(f'{self.structure_path}/relapse')
            self.relapse_ss.write(outdir=f'{self.structure_path}/relapse')

        return ctv

    def reg_structures(self, ctv):
        """ Applies the registration transform parameters to the relapse structure sets

        Parameters
        ----------
        - ctv : ROI object
            the ROI object for the planning set ctv

        Returns
        -------
        - relapse_ss_transform: StructureSet object representing the transformed relapse structure set.
        - relapse_tv_transform_roi: ROI object representing the transformed recurrence tumour volume.
        - registration_comparison_struct_metrics: DataFrame or dictionary containing registration comparison metrics.
        """

        # load patients registration parameters from output stored in results directory
        reg = Registration(f'{self.results_dir}', capture_output=True)

        # apply registration parameters to relapse ss
        relapse_ss_transform = reg.transform_structure_set(self.relapse_ss)
        relapse_ss_transform.set_image(self.adjusted_fixed)

        # write the transformed relapse structures to separate directory
        if not os.path.exists(f'{self.structure_path}/transformed_relapse'):
            relapse_ss_transform.write(
                outdir=f'{self.structure_path}/transformed_relapse')

        # get transformed recurrence tumour volume
        relapse_tv_transform_roi = relapse_ss_transform.get_roi('tumour_volume')
        relapse_tv_transform_roi.set_image(self.adjusted_fixed)

        # write dose to the results directory
        self.dose.write(f"{self.results_dir}/dose.nii")

        registration_comparison_struct_metrics = self.comparison_ss.get_comparison(
            relapse_ss_transform)

        # final view of scans with structure sets and planning dose
        # get metrics for viewing primary CTV and relapse tumour volume
        # relapse_metrics = ctv.get_comparison(relapse_tv_transform_roi)

        # self.adjusted_fixed.view(
        #   rois=[ctv, relapse_tv_transform_roi], dose=self.dose, legend=True)

        return relapse_ss_transform, relapse_tv_transform_roi, registration_comparison_struct_metrics

    def compare_structures(self, relapse_ss_transform):
        """ Compares structures between the planning structure set and the transformed relapse structure set.

        Parameters:
        ----------
        - relapse_ss_transform: StructureSet object representing the transformed relapse structure set.

        Returns:
        --------
        - comparison_metrics: DataFrame or dictionary containing comparison metrics between structures.
        """

        comparison_metrics = self.planning_ss.get_comparison(relapse_ss_transform)

        return comparison_metrics


# Create Structures
class Create_Structures(Patient_Info):
    """ A class for creating synthetic structures for a patient, inheriting from the Patient_Info class.

    Parameters:
    - p: Patient object containing relevant information.
    - trial: String representing the trial type ('import' or other).
    - recurrence_dict:  Dictionary containing naming scheme for recurrence volume.
    - ctv_names:  Dictionary containing naming scheme for Clinical Target Volumes (CTVs).
    - comparison_struct_dict:  Dictionary containing naming scheme for comparison structure information.

    Methods:
    - create_structure_spheroid(): Creates a synthetic spherical structure for a given ROI.
    - create_isodose_structures(): Creates isodose structures based on trial arm dose prescriptions.
    """

    def __init__(self, p, trial, recurrence_dict, ctv_names, comparison_struct_dict):
        """Initializes Create_Structures object with relevant parameters."""
        super().__init__(p, trial, recurrence_dict, ctv_names, comparison_struct_dict)

    def create_structure_spheroid(self, roi, im, radius=5, name=''):
        """Creates a synthetic spherical structure for a given ROI.

        Parameters:
        -----------
        - roi: ROI object representing a chosen structure.
        - im: Image object representing the image to which the structure will be added.
        - radius: Radius of the sphere that represent the registration uncertainty.
        - name: Name for the created structure.

        Returns:
        --------
        - sim_roi_sphere: ROI object representing the synthetic spherical structure.
        """
        # set image
        if im == 'adjusted_fixed':
            im = self.adjusted_fixed
        else:
            im = self.adjusted_moving

        # create synthetic image using im parameters and create synthetic sphere
        sim = SyntheticImage(shape=im.get_data().shape,
                             voxel_size=im.get_voxel_size(), origin=im.get_origin())
        sim.add_sphere(radius=radius, centre=roi.get_centroid(
            units='mm'), is_roi=True, name=name)
        sim_roi_sphere = sim.get_roi(name)

        return sim_roi_sphere

    def create_isodose_structures(self, trial_arm_dose_prescriptions):
        """ Creates isodose structures based on trial arm dose prescriptions.

        Parameters:
        -----------
        - trial_arm_dose_prescriptions: List of dose prescriptions to high, intermediate and low risk regions representing dose     
                                        prescriptions for each arm of the trial or treatment pathway

        Returns:
        ---------
        - ss: StructureSet object containing the created isodose structures.
        - trial_arm_dose_prescriptions_95: List of lists representing 95% dose prescriptions for each region for each trial arm/
                                            treatment pathway.
        """

        # get 95% of the dose prescription to each risk region
        self.trial_arm_dose_prescriptions_95 = [
            [(j*0.95) for j in i] for i in trial_arm_dose_prescriptions]

        # set colours and names for each region for visualisations
        self.color = ['blue', 'gold', 'crimson']
        self.struct = ['low_risk', 'intermediate_risk', 'high_risk']

        # create structure set object to store each of the generated isodose structures
        ss = StructureSet()
        for i in range(3):  # i = trial arm
            for j in range(3):  # j = which structure and the associated plotting colour of structure
                if self.trial_arm_dose_prescriptions_95[i][j] != 0.0:
                    roi = ROI(self.dose.get_intensity_mask(
                        self.trial_arm_dose_prescriptions_95[i][j]), color=self.color[j], name=f"isodose_{self.struct[j]}_a{i}")
                    ss.add_roi(roi)

        # set an associated image with the created 95% structure set
        ss.set_image(self.adjusted_fixed)

        return ss, self.trial_arm_dose_prescriptions_95


class Classification(Transform_Structures):
    """ A class for performing classification based on dose and volume criteria, inheriting from the Transform_Structures class.

    Parameters:
    - p: Patient object containing relevant information.
    - trial: String representing the trial type ('import' or other).
    - recurrence_dict:  Dictionary containing naming scheme for recurrence volume.
    - ctv_names:  Dictionary containing naming scheme for Clinical Target Volumes (CTVs).
    - comparison_struct_dict:  Dictionary containing naming scheme for comparison structure information.

    Methods:
    - get_centroid_dose(): obtains the dose at the centroid of a given structure.
    - get_trial_arm(): Determines the trial arm based on CTV mean dose and trial arm dose prescriptions.
    - voxtox_centroid_dose_classification(): Classifies structures based on dose at centroid of the synthetic sphere.
    - voxtox_volume_dose_classification(): Classifies structures based on how much of therecurrence sphere's volume recieves
                                        95% dose values. 
    """

    def __init__(self, p, results_dir, trial, recurrence_dict, ctv_names, comparison_struct_dict):
        """Initializes Classification object with relevant parameters."""
        super().__init__(p, results_dir, trial, recurrence_dict, ctv_names, comparison_struct_dict)

    def get_centroid_dose(self, struct):
        """ Gets the dose at the centroid of a given structure.

        Parameters:
        - struct: Structure object for which the dose is calculated.

        Returns:
        - dose_at_centroid: Dose at the structure centroid.
        """
        # get centroid of structure in voxels and obtain dose values for centroid position
        struct_centroid = list(struct.get_centroid(units='voxels').astype(int))
        dose_data = self.dose.get_data()
        dose_at_centroid = dose_data[struct_centroid[1]
                                     ][struct_centroid[0]][struct_centroid[2]]
        print(f'> Dose at {struct.name} centroid:', dose_at_centroid)

        return dose_at_centroid

    def get_trial_arm(self, trial_arm_dose_prescriptions, ctv):
        """Determines the trial arm based on CTV mean dose and trial arm dose prescriptions.

        Parameters:
        -----------
        - trial_arm_dose_prescriptions: List of lists representing dose prescriptions for each trial arm/treatment pathway.
        - ctv: ROI object representing the Clinical Target Volume.

        Returns:
        --------
        - trial_arm: Index of the determined trial arm.
        """
        ctv_mean_dose = self.dose.get_mean_dose(ctv)
        print('> CTV Mean Dose:', ctv_mean_dose)

        print(trial_arm_dose_prescriptions[2][2], trial_arm_dose_prescriptions[1]
              [2], trial_arm_dose_prescriptions[0][2])

        # create list of trial arm dose prescriptions using their index position
        # sort new list based on the ascending values
        # take the index as the trial arm, and the second value as the dose
        dose_list = [[0, trial_arm_dose_prescriptions[0][2]], [
            1, trial_arm_dose_prescriptions[1][2]], [2, trial_arm_dose_prescriptions[2][2]]]
        sorted_dose_list = sorted(dose_list, key=lambda x: x[1])
        if ctv_mean_dose >= (sorted_dose_list[2][1] - 0.5):
            trial_arm = sorted_dose_list[2][0]
        elif ctv_mean_dose >= (sorted_dose_list[1][1] - 0.5):
            trial_arm = sorted_dose_list[1][0]
        elif ctv_mean_dose >= (sorted_dose_list[0][1] - 0.5):
            trial_arm = sorted_dose_list[0][0]
        elif ctv_mean_dose <= (sorted_dose_list[0][1] - 3.0) and ctv_mean_dose >= (sorted_dose_list[2][1] + 3.0):
            trial_arm = 'check trial arms'
        else:
            trial_arm = 'nan'

        return trial_arm

    def voxtox_centroid_dose_classification(self, sim_roi_sphere, isodose_ss, trial_arm, trial_arm_dose_prescriptions):
        """ Classifies structures based on dose at the centroid of the synthetic sphere.
            If the dose at the centroid is >= the 95% dose value then it's inside the isodose structure and set == 1,
            otherwise outside the isodose structure and set == 0. 

        Parameters:
        ----------
        - sim_roi_sphere: ROI object representing the synthetic sphere.
        - isodose_ss: StructureSet object representing isodose structures.
        - trial_arm: Index of the determined trial arm.
        - trial_arm_dose_prescriptions: List of lists representing dose prescriptions for trial arm/treatment pathway.

        Returns:
        --------
        - dose_dict: Dictionary containing mean and max doses for isodose structures.
        - sim_sphere_centroid_dose: Dose at the synthetic sphere centroid.
        - centroid_dose_class_dict: Dictionary containing classification results based on dose at the centroid.
        """

        # get dose at centroid
        sim_sphere_centroid_dose = self.get_centroid_dose(sim_roi_sphere)
        # print('> Dose at rGTV centroid:', sim_sphere_centroid_dose)

        # set up empty dictionarys to store dose values
        centroid_dose_class_dict = {}
        dose_dict = {}

        # calculate mean and max dose for each isodose structure
        for isodose_roi in isodose_ss:
            isodose_mean_dose = self.dose.get_mean_dose(isodose_roi)
            isodose_max_dose = self.dose.get_max_dose_in_rois(rois=[isodose_roi])
            dose_dict[f'{isodose_roi.name}_mean'] = isodose_mean_dose
            dose_dict[f'{isodose_roi.name}_max'] = isodose_max_dose

        for isodose_roi in isodose_ss:
            if f'a{trial_arm}' in isodose_roi.name:
                # has WB and PB and TB
                if 'high_risk' in isodose_roi.name:
                    if sim_sphere_centroid_dose >= trial_arm_dose_prescriptions[trial_arm][2]:
                        # inside the structure
                        centroid_dose_class_dict[f'{isodose_roi.name}_mean'] = 1
                    else:
                        # outside the structure
                        centroid_dose_class_dict[f'{isodose_roi.name}_mean'] = 0
                elif 'intermediate_risk' in isodose_roi.name:
                    if sim_sphere_centroid_dose >= trial_arm_dose_prescriptions[trial_arm][1]:
                        # inside the structure
                        centroid_dose_class_dict[f'{isodose_roi.name}_mean'] = 1
                    else:
                        # outside the structure
                        centroid_dose_class_dict[f'{isodose_roi.name}_mean'] = 0
                else:
                    if sim_sphere_centroid_dose >= trial_arm_dose_prescriptions[trial_arm][0]:
                        # inside the structure
                        centroid_dose_class_dict[f'{isodose_roi.name}_mean'] = 1
                    else:
                        # outside the structure
                        centroid_dose_class_dict[f'{isodose_roi.name}_mean'] = 0

        return dose_dict, sim_sphere_centroid_dose, centroid_dose_class_dict

    def voxtox_volume_dose_classification(self, sim_roi_sphere, dose_dict, trial_arm, trial_arm_dose_prescriptions):
        """ Classifies structures based on how much of the recurrence sphere's volume receives 95% dose values. 

        Parameters:
        - sim_roi_sphere: ROI object representing the synthetic sphere.
        - dose_dict: Dictionary containing mean and max doses for isodose structures.
        - trial_arm: Index of the determined trial arm. 
        - trial_arm_dose_prescriptions: List of lists representing dose prescriptions for trial arm/treatment pathway.

        Returns:
        - volume_class_dict: Dictionary containing classification results based on volume and dose criteria.
        """

        # set up empty dictionary for storing volume classifications
        volume_class_dict = {}

        # determine dose mask for the recurrence sphere
        dose_in_relapse = sim_roi_sphere.get_mask().astype(int)*self.dose.get_data()
        total_voxels_in_sphere = sim_roi_sphere.get_mask().astype(int).sum()
        # print('> Total number of voxels in sphere:', total_voxels_in_sphere)

        # obtain the dose prescriptions dependent on the trial arm
        print('> Trial Arm:', trial_arm)
        if trial_arm == 0:
            isodose_95 = trial_arm_dose_prescriptions[0]
        elif trial_arm == 1:
            isodose_95 = trial_arm_dose_prescriptions[1]
        else:
            isodose_95 = trial_arm_dose_prescriptions[2]

        # for each risk region determine the number of voxels >= the 95% dose value for that risk region
        for key in dose_dict.keys():
            # print('> Key name:', key)
            if f'a{trial_arm}' in key and 'mean' in key:
                if 'high_risk' in key:
                    # print('> Structure:', key)
                    voxel_number = np.count_nonzero(
                        (dose_in_relapse >= (isodose_95[2])).astype(int))
                elif 'intermediate_risk' in key:
                    # print('> Structure:', key)
                    voxel_number = np.count_nonzero(
                        (dose_in_relapse >= (isodose_95[1])).astype(int))
                elif 'low_risk' in key:
                    # print('> Structure:', key)
                    if trial_arm == 0:
                        voxel_number = np.count_nonzero(
                            (dose_in_relapse >= (isodose_95[0])).astype(int))
                    else:
                        voxel_number = np.count_nonzero(
                            (dose_in_relapse >= (isodose_95[0])).astype(int))
                else:
                    voxel_number = 'nan'
                # print('> Voxel number:', voxel_number)

                percentage_greater = (voxel_number/total_voxels_in_sphere)*100
                # print('> Percentage greater:', percentage_greater)

                # dependent on percentage recieving the dose to each each risk region
                # allocate if sphere lies inside, outside or on the peripheral of each risk region
                if percentage_greater >= 95:
                    volume_class_dict[f'{key}'] = 1  # inside the structure

                elif percentage_greater < 95 and percentage_greater > 0:
                    # on the peripheral of the structure
                    volume_class_dict[f'{key}'] = 2

                else:
                    volume_class_dict[f'{key}'] = 0  # not in the structure

        return volume_class_dict
