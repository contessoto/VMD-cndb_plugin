/***************************************************************************
 *cr
 *cr            (C) Copyright 2020-2022 Rice University
 *cr        The Center for Theoretical Biological Physics (CTBP)
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $Author: Vinícius G. Contessoto$
 *      $Revision: 1.0 $       $Date: 2022/09/27 19:03:00 $
 *
 ***************************************************************************/

/***************************************************************************
 * 
 *      This is a modified version of H5MD VMD plugin
 *
 ***************************************************************************/

#ifndef __cndb_H
#define __cndb_H

struct cndb_file;

// all functions return 0 upon success and != 0 upon failure

//hide hdf5 error messages
void cndb_hide_hdf5_error_messages();

//show hdf5 error messages
void cndb_show_hdf5_error_messages();

/*read operations*/

// opens the file, iff it exists and creates the internal structure and goes to the first timestep
int cndb_open(struct cndb_file** _file, const char *filename, int can_open);

// close the file and frees the internal structure
int cndb_close(struct cndb_file* file);

// return the current error message 
const char* cndb_error(struct cndb_file* file);

// go to the i'th timestep
int cndb_seek_timestep(struct cndb_file* file, int i);

// reads the next timestep of all groups
int cndb_get_timestep(struct cndb_file* file, float *coords);

//reads all information regarding a given timeindependent property TODO generalize to timedependent properties, should be easy
int cndb_get_all_infromation_about_property(struct cndb_file *file, char* property, void** infromation_out);

//unfold the positions using image data
int cndb_unfold_positions(struct cndb_file* file, float* unsorted_folded_pos);

//sort data according to id datasets
int cndb_sort_data_according_to_id_datasets(struct cndb_file* file, float* to_be_sorted_data);

//reads all box informations of all groups, returns only the box information of the first group (since VMD doesn't support more than one box per file), the returned array has length 6 (3 lengths: A, B, C, and 3 angles alpha, beta, gamma)
int cndb_get_box_information(struct cndb_file* file, float* out_box_information);

//reads the number of atoms
int cndb_get_natoms(struct cndb_file* file, int* natoms);

//set the number of atoms for a dataset
int cndb_set_natoms(struct cndb_file* file, int natoms);

//reads the current time
int cndb_get_current_time(struct cndb_file* file, int* current_time);

//reads the total number of timesteps
int cndb_get_ntime(struct cndb_file* file, int* ntime);

//read time-independent dataset automatically
int cndb_read_timeindependent_dataset_automatically(struct cndb_file* file, char* dataset_name, void** _data_out, H5T_class_t* type_class_out);

//read time-independent integer dataset automatically
int cndb_read_timeindependent_dataset_int(struct cndb_file* file, char* dataset_name, int** _data_out);

//free time-independent dataset automatically
int cndb_free_timeindependent_dataset_automatically(H5T_class_t type_class, void* old_data_out, int length_array_of_strings);

//get length of one-dimensional dataset
int cndb_get_length_of_one_dimensional_dataset(struct cndb_file *file,char *dataset_name, int *length_of_dataset);

/*write operations*/

//creates a h5 file and sets the author (if not overwritten later, the username of the currently logged in user is used) and saves the reference to this bare file in a cndb_file struct, fails if a file with the provided filename already exists
int cndb_create_file(struct cndb_file **_file, const char* filename);

//sets the author's name and email. if name==NULL then the username of the currently logged in user is used, if email_address==NULL, then it is remarked, that no email was provided
int cndb_set_author(struct cndb_file* file, char* name, char* email_address);

//sets the creator program name and version. if name==NULL then the string "libcndb" is set. If version==NULL, then it is remarked, that no version was provided
int cndb_set_creator(struct cndb_file* file, char* name, char* version);

//deletes a file
int cndb_delete_file(char* filename);

//writes a dataset of given datatype
int cndb_write_dataset(struct cndb_file *file, char* absolute_name_of_dataset, hid_t datatype, void* data_in, int rank_in, hsize_t* dims_in);

//appends to an existing dataset or creates it
int cndb_append_dataset(struct cndb_file *file, char* absolute_name_of_dataset, hid_t datatype, void* positions_in, int rank_in, hsize_t* dims_in);

//gets the fill value of a dataset
int get_fill_value(struct cndb_file* file, char* absolute_name_of_dataset, void* filler);

//deletes a dataset
int cndb_delete_dataset(struct cndb_file* file, char* absolute_name_of_dataset);

#endif
