#define MAX_VAR_NAME_LENGTH     32

int exoid, exoid2, num_dim, num_nodes, num_elem, num_elem_blk;
  int *num_elem_in_block, num_node_sets, *num_nodes_per_elem, *num_attr;
  int num_side_sets, error;
  int i, j, k, m;
  int *elem_map, **connect, **ns_node_list, **ss_node_list, **node_ctr_list, **elem_list, **side_list;
  int *ebids;
  int *ssids;
  int *nsids;
  int *num_nodes_per_set, *num_elem_per_set;
  int *num_df_per_set;

  int *num_nodes_in_set, *num_elem_in_set;
  int *num_sides_in_set, *num_df_ns_in_set,*num_df_ss_in_set;
  int num_qa_rec, num_info;
  int CPU_word_size,IO_word_size;
  int num_eb_props, eb_prop_value, *eb_prop_values;
  int num_ss_props, ss_prop_value, *ss_prop_values;
  int num_ns_props, ns_prop_value, *ns_prop_values;

int num_nod_vars;
int num_time_steps;
double * time_values;

double ***nod_var_values;

  double *xcoor, *ycoor, *zcoor;
  double **ss_dist_fact;
  double **ns_dist_fact;
  
  float version, fdum;

  char *coord_names[3], *qa_record[100][4], **info, **nodal_var_names;
  char title[MAX_LINE_LENGTH+1], elem_type[100][MAX_STR_LENGTH+1];
  char *eb_prop_names[3],*ss_prop_names[3],*ns_prop_names[3];
  char *cdum;
