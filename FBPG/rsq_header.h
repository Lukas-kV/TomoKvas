typedef struct s_ima_data_type {
/*---------------------------------------------*/
char check[16];
int data_type;
int nr_of_bytes; /* either one of them */
int nr_of_blocks; /* or both, but min. of 1 */
int patient_index; /* 1 block = 512 bytes */
int scanner_id;
int creation_date[2];
/*---------------------------------------------*/
int dimx_p;
int dimy_p;
int dimz_p;
int dimx_um;
int dimy_um;
int dimz_um;
int slice_thickness_um;
int slice_increment_um;
int slice_1_pos_um;
int min_data_value;
int max_data_value;
int mu_scaling; /* p(x,y,z)/mu_scaling = value [1/cm] */
int nr_of_samples;
int nr_of_projections;
int scandist_um;
int scanner_type;
int sampletime_us;
int index_measurement;
int site; /* Coded value */
int reference_line_um;
int recon_alg; /* Coded value */
char name[40];
int energy; /*V */
int intensity; /* uA */
int fill[83];
/*---------------------------------------------*/
int data_offset; /* in 512-byte-blocks */
} ima_data_type, *ima_data_typeP; 
