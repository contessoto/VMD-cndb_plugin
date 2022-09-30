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
 *      All credits to the authors named below
 *
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $Author: Jonas Landsgesell, Sascha Ehrhardt S$
 *      $Revision: 1.5 $       $Date: 2014/03/27 19:03:00 $
 *
 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "molfile_plugin.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "libcndb.h"

//element symbols
static const char *element_symbols[] = { 
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P" , "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", 
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg"
};

const char default_name[16]="C";
const char default_type[16]="C";
const char default_resname[8]="";
const int default_resid= 1;
const char default_types[16]="NA";
const char default_segid[8]= "";
const char default_chain[2]= "";
const char default_altloc[2]= "";
const char default_insertion[2]= "";
const float default_occupancy= 1.0;
const float default_bfactor= 1.0; 
const float default_mass= 1.0;
const float default_charge= 0.0;
const float default_radius= 0.5;
const int default_atomicnumber= 1;

static void *open_cndb_read(const char *filename, const char *filetype, int *natoms){
	cndb_hide_hdf5_error_messages();
	struct cndb_file* file;
	cndb_open(&file,filename,-1);
	cndb_get_natoms(file, natoms);
	return file;
}


/* read the coordinates */
static int read_cndb_timestep(void *_file, int natoms, molfile_timestep_t *ts) {
	struct cndb_file* file=_file;
	int status=MOLFILE_SUCCESS;
	int current_time;
	cndb_get_current_time(file,&current_time);
	int ntime;
	cndb_get_ntime(file,&ntime);
	if(current_time>=ntime)
		return MOLFILE_ERROR;
	if (ts != NULL) { //skip reading if ts is NULL pointer (needs modification of the timestep below)
		cndb_get_natoms(file, &natoms);

		//read boxinformation
		float box_information[6];
		int status_box=cndb_get_box_information(file, box_information);
		if(status_box==0){
			ts->A=box_information[0];
			ts->B=box_information[1];
			ts->C=box_information[2];
			ts->alpha=box_information[3];
			ts->beta=box_information[4];
			ts->gamma=box_information[5];
		}
		ts->velocities = NULL;
		ts->physical_time = 0.0;
		//read coords
		cndb_get_timestep(file, ts->coords);
		cndb_unfold_positions(file, ts->coords);
		cndb_sort_data_according_to_id_datasets(file, ts->coords);

	}
	
	int status_seek=cndb_seek_timestep(file, current_time+1); //modify timestep in the internal state of the plugin for this file
	if(status_seek!=0){
		status= MOLFILE_SUCCESS;
	}else if(status_seek!=0 ){
		status= MOLFILE_ERROR;
	}
	return status;
}

//find index of species "int species" in data_index_species which has length len_index
int find_index_of_species(int *data_index_species, int species,int len_index){
	int ndx=-1;
	for(int i=0;i<len_index;i++){
		if(species== data_index_species[i]){
			ndx=i;
			break;
		}
	}

	return ndx;
}

//count the different species that are present in data_species, requires pointer to array and length of array as argument
int find_first_index_in_array(int *data_species, int index_of_current_value, int length){
	int index_at_which_the_value_of_current_index_first_occured=-1;
	for(int i=0;i<length;i++){
		if(data_species[i]==data_species[index_of_current_value]){
			index_at_which_the_value_of_current_index_first_occured=i;
			break;
		}
	}
	return index_at_which_the_value_of_current_index_first_occured;	
}

int count_different_species(int *data_species, int length){
	int count =0;
	for(int i=0;i<length;i++){
		if(find_first_index_in_array(data_species,i,length)==i){
			count+=1; //add plus 1 if the first occurance of a species happens at the current index i since then a new species is found
		}	
	}
	return count;
}



//check for consistency between the number of species and the corresponding index of species
int check_consistency_species_index_of_species(struct cndb_file *file, int len_data_index_species, int* data_species){
	int natoms;
	cndb_get_natoms(file, &natoms);
	if(count_different_species(data_species,natoms) == len_data_index_species){
		return 0;
	}else{
		return -1;
	}
}


//load whole VMD structure
int read_cndb_structure_vmd_structure(void *_file, int *optflags,molfile_atom_t *atoms) {
	molfile_atom_t *atom;
	*optflags = MOLFILE_ATOMICNUMBER | MOLFILE_MASS | MOLFILE_RADIUS |  MOLFILE_INSERTION | MOLFILE_CHARGE; // we read the optional attributes atomicnumber, mass, radius and charge
	struct cndb_file* file=_file;
	int natoms;
	cndb_get_natoms(file, &natoms);

	printf("Natoms %d \n", natoms);

	//load index of species
	int* data_index_species;
	H5T_class_t type_class_index_species;
	int status_index_species=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/indexOfSpecies", (void**) &data_index_species, &type_class_index_species);

	//load species
	float* data_species_float;
	H5T_class_t type_class_species=H5T_FLOAT;
	char* species_property="species";
	int status_read_species=cndb_get_all_infromation_about_property(file,species_property ,(void**) &data_species_float);
	int data_species[natoms];
	if(status_read_species>=0){
		for(int i=0;i<natoms;i++){
			data_species[i]=(int) data_species_float[i];	
		}
	}
	//load mass
	float* data_mass;
	H5T_class_t type_class_mass=H5T_FLOAT;
	char* mass_property="mass";
	int status_read_mass=cndb_get_all_infromation_about_property(file, mass_property, (void**) &data_mass);


	//load charge
	float* data_charge;
	H5T_class_t type_class_charge=H5T_FLOAT;
	char* charge_property="charge";
	int status_read_charge=cndb_get_all_infromation_about_property(file,charge_property , (void**) &data_charge);

	//Declaring variables here since if one would declare them later in the else branch one could not access them to free them later, after the atoms have been assigned to their values
	char **data_name;
	int status_read_name=-1;
	H5T_class_t type_class_name;
	char **data_type;
	int status_read_type=-1;
	H5T_class_t type_class_type;
	float* data_radius;
	H5T_class_t type_class_radius;
	int status_read_radius=-1;
	int* data_atomicnumber;
	H5T_class_t type_class_atomicnumber;
	int status_read_atomicnumber=-1;	
	char** data_segid;
	int status_read_segid=-1;
	H5T_class_t type_class_segid;
	int* data_resid;
	int status_read_resid=-1;
	H5T_class_t type_class_resid;
	char** data_resname;
	int status_read_resname=-1;
	H5T_class_t type_class_resname;
	char** data_chain;
	int status_read_chain=-1; 
	H5T_class_t type_class_chain;

	char **data_types;
	int status_read_types=-1;
	H5T_class_t type_class_types;
	
	int len_data_index_species;
	int len_data_resname;
	int species_check = -1;
	if(status_read_species>=0){
		cndb_get_length_of_one_dimensional_dataset(file,"/parameters/vmd_structure/indexOfSpecies",&len_data_index_species);
		species_check=check_consistency_species_index_of_species(file, len_data_index_species, data_species);
	}
	if(species_check!=0){
		// printf("NOTE: /parameters/vmd_structure/indexOfSpecies does not contain as much different species as there are species present in the different groups /particles/*/species !\n");
		// printf("Skipping index of species related data.\n");
		// status_read_types=cndb_read_timeindependent_dataset_automatically(file, "types",(void**) &data_types, &type_class_types);
		// printf("status_read_types - if %d \n", status_read_types);
		printf("...\n ...\n");
	}else{
		status_read_name=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/name",(void**) &data_name, &type_class_name);
		printf("status_read_name  %d \n", status_read_name);
		status_read_type=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/type",(void**) &data_type, &type_class_type);
		status_read_radius=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/radius",(void**) &data_radius, &type_class_radius);
		status_read_atomicnumber=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/atomicnumber",(void**) &data_atomicnumber, &type_class_atomicnumber);
	}
	// Added to read chromatin types folder
	status_read_types=cndb_read_timeindependent_dataset_automatically(file, "types",(void**) &data_types, &type_class_types);
	printf("status_read_types - out else %d \n", status_read_types);
	status_read_segid=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/segid",(void**) &data_segid, &type_class_segid);
	//load resid
	printf("status_read_segid %d \n", status_read_segid);
	// status_read_resid=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/resid",(void**) &data_resid, &type_class_resid);
	status_read_resid=0;
	printf("status_read_resid %d \n", status_read_resid);
	if(status_read_resid==0){
		status_read_resname=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/resname",(void**) &data_resname, &type_class_resname);
		status_read_chain=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/chain",(void**) &data_chain, &type_class_chain);
		cndb_get_length_of_one_dimensional_dataset(file,"types",&len_data_resname);
		printf("len_data_resname %d \n", len_data_resname);
	}
	//give data to VMD
	for (int i = 0; i < natoms; i++) {
		atom = atoms + i;
		// printf("atom %s, i %d \n", atoms,i);
		//set species related properties
		int index_of_species=-1;
		if(status_index_species==0 && status_read_species>=0)
			index_of_species=find_index_of_species(data_index_species,data_species[i],len_data_index_species);
		if(status_read_type==0 && status_index_species==0 &&index_of_species>=0)
			strncpy(atom->type, data_type[index_of_species], 16*sizeof(char));	//set type for atom of species
			// strncpy(atom->type, default_type, 16*sizeof(char));	//set type for atom of species
		else
			strncpy(atom->type,default_type,16*sizeof(char));
		if(status_read_atomicnumber==0 && status_index_species==0 &&index_of_species>=0){	//set atomicnumber
			if(data_atomicnumber[index_of_species]<112 ){
					atom->atomicnumber = data_atomicnumber[index_of_species]; 	
			}else{
				atom->atomicnumber = data_atomicnumber[index_of_species]%112;
				printf("ERROR: Too big atomic number was set in dataset. Using modulo %%112 to produce smaller atomic number");
				return MOLFILE_ERROR;
			}
		}else{
			if(status_read_species>=0){
				atom->atomicnumber = data_species[i]%112;
			}else{
				atom->atomicnumber = default_atomicnumber;
			}
		}

		// printf("atom %d, i %d, atom->type %s, default_type %s \n", atom->atomicnumber,i,atom->type,default_type);



		if (status_read_name==0 && status_index_species==0 && index_of_species>=0){
			strncpy(atom->name,data_name[index_of_species],16*sizeof(char));	//set elementname for atom of species
		}
		else{
			strncpy(atom->name,element_symbols[atom->atomicnumber],16*sizeof(char));
		}

		// printf("atom %s, i %d, atom->type %s, element_symbols %s \n", atom->name,i,atom->type,element_symbols[atom->atomicnumber]);

		if(status_read_mass==0 && status_index_species==0 && index_of_species>=0)
			atom->mass = data_mass[i];	//set mass for atom of with id i
		else
			atom->mass=default_mass;
		if(status_read_radius==0 && status_index_species==0 && index_of_species>=0)	
			atom->radius = data_radius[index_of_species];	//set radius for atom of species
		else
			atom->radius=default_radius;
		if(status_read_charge==0)
			atom->charge=data_charge[i];	//set charge
		else
			atom->charge=default_charge;
		if(status_read_segid==0)
			strncpy(atom->segid,data_segid[i],sizeof(8));	//set segid (segments may also consist of pure solvent or lipid)
		else
			strncpy(atom->segid,default_segid,sizeof(8));
		if(status_read_resid==0)
			atom->resid=data_resid[i];	//set resid
		else
			atom->resid=default_resid;
		if(status_read_resname==0)		
			strncpy(atom->resname,data_resname[data_resid[i]],8*sizeof(char));	//set resname
		else
			strncpy(atom->resname,default_resname,8*sizeof(char));
		if(status_read_chain==0)
			strncpy(atom->chain,data_chain[data_resid[i]],2*sizeof(char));	//set chain
		else
			strncpy(atom->chain,default_chain,2*sizeof(char));

		if(status_read_types==0)
			strncpy(atom->resname,data_types[i],2*sizeof(char));	//set Types
		else
			strncpy(atom->resname,default_types,2*sizeof(char));

		// printf("atom %d, i %d, atom->resname %s, element_symbols - Here %s \n", atom->resid,i,atom->resname,atom->chain);
		printf("atom %d, i %d, atom->resname %s, default_types %s - data_chain[data_resid[i] %d\n", atom->resid,i,atom->resname,default_types,status_read_types);
	}

	


	//After assignment free used resources
	if(status_index_species==0)
		cndb_free_timeindependent_dataset_automatically(type_class_index_species,data_index_species, 0);
	if(status_read_species>=0)	
		cndb_free_timeindependent_dataset_automatically(type_class_species,data_species_float, 0);
	if(status_read_name==0)
		cndb_free_timeindependent_dataset_automatically(type_class_name,data_name, len_data_index_species);
	if(status_read_type==0)
		cndb_free_timeindependent_dataset_automatically(type_class_type,data_type, len_data_index_species);
	if(status_read_mass==0)
		cndb_free_timeindependent_dataset_automatically(type_class_mass,data_mass, 0);
	if(status_read_radius==0)
		cndb_free_timeindependent_dataset_automatically(type_class_radius,data_radius, 0);
	if(status_read_atomicnumber==0)
		cndb_free_timeindependent_dataset_automatically(type_class_atomicnumber,data_atomicnumber, 0);
	if(status_read_charge==0)
		cndb_free_timeindependent_dataset_automatically(type_class_charge,data_charge, 0);
	if(status_read_segid==0)
		cndb_free_timeindependent_dataset_automatically(type_class_segid,data_segid, natoms);
	if(status_read_resid==0)
		cndb_free_timeindependent_dataset_automatically(type_class_resid,data_resid, 0);
	if(status_read_resname==0)
		cndb_free_timeindependent_dataset_automatically(type_class_resname,data_resname, len_data_resname);
	if(status_read_chain==0)
		cndb_free_timeindependent_dataset_automatically(type_class_chain,data_chain, len_data_resname);


	return MOLFILE_SUCCESS;
}


static int cndb_get_bonds(void *_file, int *nbonds, int **from, int **to, float **bondorder, int **bondtype,  int *nbondtypes, char ***bondtypename){
	*bondorder = NULL;
	*bondtype=NULL;
	*nbondtypes=0;
	bondtypename=NULL;
	struct cndb_file* file=_file; 
	H5T_class_t type_class_bond_from;
	H5T_class_t type_class_bond_to;
	
	cndb_get_length_of_one_dimensional_dataset(file,"/parameters/vmd_structure/bond_from",nbonds); //save number of bonds to *nbonds
	
	
	int status_read_bond_from=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/bond_from",(void**) from, &type_class_bond_from);
	int status_read_bond_to=cndb_read_timeindependent_dataset_automatically(file, "/parameters/vmd_structure/bond_to",(void**) to, &type_class_bond_to);
	
	if(status_read_bond_from==0 && status_read_bond_to ==0){ 
		printf("read %d bonds \n",*nbonds);
  		return MOLFILE_SUCCESS;
	}else{
		return MOLFILE_ERROR;		
	}
}

void close_file(void* _file){
	struct cndb_file* file=_file;
	cndb_close(file);
	//TODO bonds from, to need to be freed, compare line 01050 http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/psfplugin_8c-source.html
}


//write operations
void* open_cndb_write(const char* filename, const char *filetype, int natoms){
	struct cndb_file* file;
	cndb_create_file(&file, filename);
	// cndb_set_natoms(file, natoms);
	return file;
}

static int write_cndb_timestep(void* _file, const molfile_timestep_t *ts){
	struct cndb_file* file=(struct cndb_file*) _file;
	int natoms;
	cndb_get_natoms(file, &natoms);
	int rank_in=3;
	hsize_t dims_in[3];
	dims_in[0] =1;	//for one timestep
	dims_in[1] =natoms;
	dims_in[2] =3; //number of space dimensions

	int status= cndb_append_dataset(file, "/particles/atoms/position/value", H5T_NATIVE_FLOAT, (void*) ts->coords, rank_in, dims_in);
	
	if(status==0)
		return MOLFILE_SUCCESS;
	else
		return MOLFILE_ERROR;
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
	memset(&plugin, 0, sizeof(molfile_plugin_t));
	plugin.abiversion = vmdplugin_ABIVERSION;
	plugin.type = MOLFILE_PLUGIN_TYPE;
	plugin.name = "cndb";
	plugin.prettyname = "cndb";
	plugin.author = "Sascha Ehrhardt, Jonas Landsgesell";
	plugin.majorv = 1;
	plugin.minorv = 5;
	plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
	plugin.filename_extension = "cndb";
	plugin.open_file_read = open_cndb_read;
	plugin.read_structure = read_cndb_structure_vmd_structure;
	plugin.read_next_timestep = read_cndb_timestep;
	plugin.read_bonds = cndb_get_bonds;
	plugin.close_file_read = close_file;
	plugin.open_file_write = open_cndb_write;
	//plugin.write_structure = write_cndb_vmd_structure;
	plugin.write_timestep = write_cndb_timestep;
	plugin.close_file_write = close_file;
	

	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
	(*cb)(v, (vmdplugin_t *) &plugin);
	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
	return VMDPLUGIN_SUCCESS;
}
