/* exodus reader gets all infor from exodus file
   derived from Larry Schoof's testrdwt program  */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "exodusII.h"
#include "netcdf.h"


#include "exo_rw_data.h"

FILE *log_file;

int get_exodus (char* input_file)

{

  log_file = fopen("exo.log","w");
  if (log_file == NULL)
    {
      printf("error opening log file\n");
    }


  cdum = 0;

  /* ex_opts(2); */

  /* Specify compute and i/o word size */

  CPU_word_size = sizeof(double);			/* sizeof(float) */
  IO_word_size = 8;			/* float */

  /* open EXODUS II file for reading */

  exoid = ex_open (input_file,		/* filename path */
		   EX_READ,		/* access mode */
		   &CPU_word_size,	/* CPU float word size in bytes */
		   &IO_word_size,	/* I/O float word size in bytes */
		   &version);		/* returned version number */
  fprintf (log_file,"after ex_open for test.exo\n");
  fprintf (log_file," cpu word size: %d io word size: %d\n",CPU_word_size,IO_word_size);



  /* read initialization parameters */

  error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem,
		       &num_elem_blk, &num_node_sets, &num_side_sets);

  fprintf (log_file,"after ex_get_init, error = %d\n", error);
  

  /* read nodal coordinate values */

  xcoor = (double *) calloc(num_nodes, sizeof(double));
  ycoor = (double *) calloc(num_nodes, sizeof(double));
  if (num_dim >= 3)
    zcoor = (double *) calloc(num_nodes, sizeof(double));
  else
    zcoor = 0;
 
  error = ex_get_coord (exoid, xcoor, ycoor, zcoor);
  fprintf (log_file,"\nafter ex_get_coord, error = %3d\n", error);
 

 
  /* read nodal coordinate names */

  for (i=0; i<num_dim; i++)
    {
      coord_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }
 
  error = ex_get_coord_names (exoid, coord_names);
  fprintf (log_file,"\nafter ex_get_coord_names, error = %3d\n", error);
 


  /* read element order map */

  elem_map = (int *) calloc(num_elem, sizeof(int));
 
  error = ex_get_map (exoid, elem_map);
  fprintf (log_file,"\nafter ex_get_map, error = %3d\n", error);
 


  /* read element block parameters and element connectivity */

  ebids = (int *) calloc(num_elem_blk, sizeof(int));
  error = ex_get_elem_blk_ids (exoid, ebids);
  fprintf (log_file,"\nafter ex_get_elem_blk_ids, error = %3d\n", error);

  num_elem_in_block    = (int *) calloc(num_elem_blk, sizeof(int));
  num_nodes_per_elem   = (int *) calloc(num_elem_blk, sizeof(int));
  num_attr             = (int *) calloc(num_elem_blk, sizeof(int));
  connect              = (int **) calloc(num_elem_blk, sizeof(int*));
  /*   *elem_type            = (char **) calloc(num_elem_blk,sizeof(char *)); */

  for (i=0; i<num_elem_blk; i++)
    {
      error = ex_get_elem_block (exoid, ebids[i], elem_type[i],
				 &num_elem_in_block[i],
				 &num_nodes_per_elem[i], &num_attr[i]);
      fprintf (log_file,"\nafter ex_get_elem_block, error = %d\n", error);
 
      connect[i] = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]),
				  sizeof(int));
 
      error = ex_get_elem_conn (exoid, ebids[i], connect[i]);
      fprintf (log_file,"\nafter ex_get_elem_conn, error = %d\n", error);
    }

  /* read block properties */

  error = ex_inquire (exoid, EX_INQ_EB_PROP, &num_eb_props, &fdum, cdum);
  fprintf (log_file,"\nafter ex_inquire, error = %d\n", error);

  for (i=0; i<num_eb_props; i++)
    {
      eb_prop_names[i] = (char *) calloc ((MAX_VAR_NAME_LENGTH+1), sizeof(char));
    }
 
  error = ex_get_prop_names(exoid,EX_ELEM_BLOCK,eb_prop_names);
  fprintf (log_file,"after ex_get_prop_names, error = %d\n", error);
 
  for (i=0; i<num_eb_props; i++)
    {
      for (j=0; j<num_elem_blk; j++)
	{
	  error = ex_get_prop(exoid, EX_ELEM_BLOCK, ebids[j], eb_prop_names[i],
			      &eb_prop_value);
	  fprintf (log_file,"after ex_get_prop, error = %d\n", error);
	}
    }

  /* read individual node sets */

  nsids              = (int *) calloc(num_node_sets, sizeof(int));
  num_nodes_in_set   = (int *) calloc(num_node_sets, sizeof(int));
  num_df_ns_in_set      = (int *) calloc(num_node_sets, sizeof(int));
  ns_node_list       = (int **) calloc(num_node_sets,sizeof(int*));
  ns_dist_fact       = (double **) calloc(num_node_sets,sizeof(double*));
  error = ex_get_node_set_ids (exoid, nsids);
  fprintf (log_file,"\nafter ex_get_node_set_ids, error = %3d\n", error);
 
  for (i=0; i< num_node_sets; i++)
    {
      error = ex_get_node_set_param (exoid, nsids[i],
				     &num_nodes_in_set[i], &num_df_ns_in_set[i]);
      fprintf (log_file,"\nafter ex_get_node_set_param, error = %3d\n", error);

      ns_node_list[i] = (int *) calloc(num_nodes_in_set[i], sizeof(int));
      ns_dist_fact[i] = (double *) calloc(num_df_ns_in_set[i] , sizeof(double));
       
      error = ex_get_node_set (exoid, nsids[i], ns_node_list[i]);

      fprintf (log_file,"\nafter ex_get_node_set, error = %3d\n", error);

      if (num_df_ns_in_set[i] > 0)
	{
	  error = ex_get_node_set_dist_fact (exoid, nsids[i], ns_dist_fact[i]);
	  fprintf (log_file,"\nafter ex_get_node_set_dist_fact, error = %3d\n", error);
	}
    }

  /* read node set properties */
  error = ex_inquire (exoid, EX_INQ_NS_PROP, &num_ns_props, &fdum, cdum);
  fprintf (log_file,"\nafter ex_inquire, error = %d\n", error);
 
  for (i=0; i<num_ns_props; i++)
    {
      ns_prop_names[i] = (char *) calloc ((MAX_VAR_NAME_LENGTH+1), sizeof(char));
    }
  ns_prop_values = (int *) calloc (num_node_sets, sizeof(int));
 
  error = ex_get_prop_names(exoid,EX_NODE_SET,ns_prop_names);
  fprintf (log_file,"after ex_get_prop_names, error = %d\n", error);

  for (i=0; i<num_ns_props; i++)
    {
      error = ex_get_prop_array(exoid, EX_NODE_SET, ns_prop_names[i],
				ns_prop_values);
      fprintf (log_file,"after ex_get_prop_array, error = %d\n", error);
    }

  /* read individual side sets */
  if(num_side_sets > 0)
    {
      ssids = (int *) calloc(num_side_sets, sizeof(int));
      num_sides_in_set = (int *) calloc(num_side_sets, sizeof(int));
      num_df_ss_in_set = (int *) calloc(num_side_sets, sizeof(int));
      num_elem_in_set = (int *) calloc(num_side_sets, sizeof(int));
      ss_dist_fact = (double **) calloc(num_side_sets,sizeof(double*));
      ss_node_list = (int **) calloc(num_side_sets,sizeof(int*));
      elem_list    = (int **) calloc(num_side_sets,sizeof(int*));
      side_list    = (int **) calloc(num_side_sets,sizeof(int*));
      node_ctr_list    = (int **) calloc(num_side_sets,sizeof(int*));
      error = ex_get_side_set_ids (exoid, ssids);
      fprintf (log_file,"\nafter ex_get_side_set_ids, error = %3d\n", error);
 
      for (i=0; i<num_side_sets; i++)
	{
	  error = ex_get_side_set_param (exoid, ssids[i], &num_sides_in_set[i],
					 &num_df_ss_in_set[i]);
	  fprintf (log_file,"\nafter ex_get_side_set_param, error = %3d\n", error);
	  
	  /* Note: The # of elements is same as # of sides!  */
	  num_elem_in_set[i] = num_sides_in_set[i];
	  elem_list[i] = (int *) calloc(num_elem_in_set[i], sizeof(int));
	  side_list[i] = (int *) calloc(num_sides_in_set[i], sizeof(int));
	  node_ctr_list[i] = (int *) calloc(num_elem_in_set[i], sizeof(int));
	  ss_node_list[i] = (int *) calloc(num_elem_in_set[i]*21, sizeof(int));
	  ss_dist_fact[i] = (double *) calloc(num_df_ss_in_set[i], sizeof(double));
 
	  error = ex_get_side_set (exoid, ssids[i], elem_list[i], side_list[i]);
	  fprintf (log_file,"\nafter ex_get_side_set, error = %3d\n", error);
 
	  error = ex_get_side_set_node_list (exoid, ssids[i], node_ctr_list[i],
					     ss_node_list[i]);
	  fprintf (log_file,"\nafter ex_get_side_set_node_list, error = %3d\n", error);
 
	  if (num_df_ss_in_set[i] > 0)
	    {
	      error = ex_get_side_set_dist_fact (exoid, ssids[i], ss_dist_fact[i]);
	      fprintf (log_file,"\nafter ex_get_side_set_dist_fact, error = %3d\n", error);
	    }
	}

      /* read side set properties */
      error = ex_inquire (exoid, EX_INQ_SS_PROP, &num_ss_props, &fdum, cdum);
      fprintf (log_file,"\nafter ex_inquire, error = %d\n", error);
 
      for (i=0; i<num_ss_props; i++)
	{
	  ss_prop_names[i] = (char *) calloc ((MAX_VAR_NAME_LENGTH+1), sizeof(char));
	}
 
      error = ex_get_prop_names(exoid,EX_SIDE_SET,ss_prop_names);
      fprintf (log_file,"after ex_get_prop_names, error = %d\n", error);
 
 
      for (i=0; i<num_ss_props; i++)
	{
	  for (j=0; j<num_side_sets; j++)
	    {
	      error = ex_get_prop(exoid, EX_SIDE_SET, ssids[j], ss_prop_names[i],
				  &ss_prop_value);
	      fprintf (log_file,"after ex_get_prop, error = %d\n", error);

	    }
	}

    }

  /* read QA records */

  ex_inquire (exoid, EX_INQ_QA, &num_qa_rec, &fdum, cdum);
  for (i=0; i < num_qa_rec; i++)
    {
      for (j=0; j<4; j++)
	{
	  qa_record[i][j] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
	}
    }

  error = ex_get_qa (exoid, qa_record);
  fprintf (log_file,"\nafter ex_get_qa, error = %3d\n", error);


  /* read information records */

  error = ex_inquire (exoid, EX_INQ_INFO, &num_info, &fdum, cdum);
  fprintf (log_file,"\nafter ex_inquire, error = %3d\n", error);

  if(num_info > 0)
    {
      info = (char **) calloc(num_info, sizeof(char*));
      for (i=0; i<num_info; i++)
	{
	  info[i] = (char *) calloc ((MAX_LINE_LENGTH+1), sizeof(char));
	}

      error = ex_get_info (exoid, info);
      fprintf (log_file,"\nafter ex_get_info, error = %3d\n", error);
    }

  
  error = ex_get_var_param (exoid, "n", &num_nod_vars);
  fprintf (log_file,"after ex_get_var_param, error = %d %d\n", error, num_nod_vars);


  nodal_var_names = (char**)calloc(num_nod_vars,sizeof(char*));
  for(i = 0; i < num_nod_vars;i++){
    nodal_var_names[i] = (char*)calloc((MAX_STR_LENGTH+1),sizeof(char));
  }

  error = ex_get_var_names (exoid, "n", num_nod_vars, nodal_var_names);
  fprintf (log_file,"after ex_get_var_names, error = %d\n", error);



 /* determine how many time steps are stored */

  error = ex_inquire (exoid, EX_INQ_TIME, &num_time_steps, &fdum, cdum);
  fprintf (log_file,"\nafter ex_inquire, error = %3d\n", error);
  fprintf (log_file,"There are %2d time steps in the database.\n", num_time_steps);
  
  
  time_values = (double *) calloc (num_time_steps, sizeof(double));
  
  error = ex_get_all_times (exoid, time_values);
  fprintf (log_file,"\nafter ex_get_all_times, error = %3d\n", error);
  

   nod_var_values = (double ***)calloc(num_time_steps,sizeof(double **));
   for(i = 0; i < num_time_steps; i++){
     nod_var_values[i] = (double **)calloc(num_nod_vars,sizeof(double *));
     for(j = 0;j < num_nod_vars;j++){
       nod_var_values[i][j] = (double *)calloc(num_nodes,sizeof(double));
       error = ex_get_nodal_var(exoid, i+1, j+1, num_nodes, nod_var_values[i][j]);
       fprintf (log_file,"\nafter ex_get_nodal_var, error = %3d\n", error);
     }
   }



  /* close the EXODUS files */

  error = ex_close (exoid);
  fprintf (log_file,"after ex_close, error = %d\n", error);
  fflush(log_file);
  return 0;
}

int put_exodus(char* output_file)
{
  if (log_file == NULL)
    {
      log_file = fopen("exo.log","w");
    }
  if (log_file == NULL)
    {
      fprintf(log_file,"error opening log file\n");
    }
  
  /* create EXODUS II file for writing */

  exoid2= ex_create (output_file,	/* filename path */
		     EX_CLOBBER,	/* create mode */
		     &CPU_word_size,	/* CPU float word size in bytes */
		     &IO_word_size);	/* I/O float word size in bytes */
  fprintf (log_file,"after ex_create for %s, exoid = %d\n",output_file, exoid2);

  /* write initialization parameters */

  error = ex_put_init (exoid2, title, num_dim, num_nodes, num_elem,
		       num_elem_blk, num_node_sets, num_side_sets);

  fprintf (log_file,"after ex_put_init, error = %d\n", error);

  /* write nodal coordinate values */

  error = ex_put_coord (exoid2, xcoor, ycoor, zcoor);
  fprintf (log_file,"after ex_put_coord, error = %d\n", error);


  /* write nodal coordinate names */

  error = ex_put_coord_names (exoid2, coord_names);
  fprintf (log_file,"after ex_put_coord_names, error = %d\n", error);


  /* write element order map */

  error = ex_put_map (exoid2, elem_map);
  fprintf (log_file,"after ex_put_map, error = %d\n", error);


  /* write element block parameters and element connectivity */
   
  for (i=0; i<num_elem_blk; i++)
    {

      error = ex_put_elem_block (exoid2, ebids[i], elem_type[i], num_elem_in_block[i],
				 num_nodes_per_elem[i], num_attr[i]);
      fprintf (log_file,"after ex_put_elem_block, error = %d\n", error);
 
      error = ex_put_elem_conn (exoid2,ebids[i], connect[i]);
      fprintf (log_file,"after ex_put_elem_conn, error = %d\n", error);

      /*       free (connect); */
    }






  /* write element block properties */
  if(num_eb_props > 1)
    {
      error = ex_put_prop_names(exoid2,EX_ELEM_BLOCK,num_eb_props,eb_prop_names);
      fprintf (log_file,"after ex_put_prop_names, error = %d\n", error);
      
      for (i=0; i<num_eb_props; i++)
	{
	  for (j=0; j<num_elem_blk; j++)
	    {
	      if (i>0) {   /* first property is the ID which is already stored */
		error = ex_put_prop(exoid2, EX_ELEM_BLOCK, ebids[j], eb_prop_names[i], 
				    eb_prop_value);
		fprintf (log_file,"after ex_put_prop, error = %d\n", error);
	      }
	    }
	}

    }






  /* write individual node sets */

  for (i=0; i < num_node_sets; i++)
    {
      error = ex_put_node_set_param (exoid2, nsids[i], num_nodes_in_set[i], 
				     num_df_ns_in_set[i]);
fprintf(log_file,"\nnsids %d, num_nodes_in_set %d, num_df_ns_in_set %d",nsids[i], num_nodes_in_set[i],num_df_ns_in_set[i]);
      fprintf (log_file,"after ex_put_node_set_param, error = %d\n", error);

      error = ex_put_node_set (exoid2, nsids[i], ns_node_list[i]);
      fprintf (log_file,"after ex_put_node_set, error = %d\n", error);

      if (num_df_ns_in_set[i] > 0)
	{
	  error = ex_put_node_set_dist_fact (exoid2, nsids[i], ns_dist_fact[i]);
	  fprintf (log_file,"after ex_put_node_set, error = %d\n", error);
	}
 
      /*       free (node_list); */
      /*       free (ns_dist_fact); */
    }



  if((num_ns_props > 0) && (num_node_sets > 0))
    {
      /* write node set properties */
      error = ex_put_prop_names(exoid2,EX_NODE_SET,num_ns_props,ns_prop_names);
      fprintf (log_file,"after ex_put_prop_names, error = %d\n", error);
 
      for (i=0; i<num_ns_props; i++)
	{
	  error = ex_put_prop_array(exoid2, EX_NODE_SET, ns_prop_names[i], ns_prop_values);
	  fprintf (log_file,"after ex_put_prop_array, error = %d\n", error);
	}   
    }


  /* write individual side sets */
  if(num_side_sets > 0)
    {
      for (i=0; i<num_side_sets; i++)
	{
	  error = ex_put_side_set_param (exoid2, ssids[i], num_sides_in_set[i], 
					 num_df_ss_in_set[i]);
	  fprintf (log_file,"after ex_put_side_set_param, error = %d\n", error);
 
	  error = ex_put_side_set (exoid2, ssids[i], elem_list[i], side_list[i]);
	  fprintf (log_file,"after ex_put_side_set, error = %d\n", error);

	  if (num_df_ss_in_set[i] > 0)
	    {
	      error = ex_put_side_set_dist_fact (exoid2, ssids[i], ss_dist_fact[i]);
	      fprintf (log_file,"after ex_put_side_set_dist_fact, error = %d\n", error);
	    }
 

	}

	  
      for (i=0; i<num_ss_props; i++)
	{
	  for (j=0; j<num_side_sets; j++)
	    {
	      if (i>0) {  /* first property is ID so it is already stored */
		error = ex_put_prop(exoid2, EX_SIDE_SET, ssids[j], ss_prop_names[i], 
				    ss_prop_value);
		fprintf (log_file,"after ex_put_prop, error = %d\n", error);
	      }
	    }
	}


    }

  if (num_qa_rec > 0)
    {
      error = ex_put_qa (exoid2, num_qa_rec, qa_record);
      fprintf (log_file,"after ex_put_qa, error = %d\n", error);
    }  
  if(num_info > 0)
    {
      error = ex_put_info (exoid2, num_info, info);
      fprintf (log_file,"after ex_put_info, error = %d\n", error);


    }

  if(num_nod_vars > 0){
    error = ex_put_var_param (exoid, "n", num_nod_vars);
    fprintf (log_file,"after ex_put_var_param, error = %d %d\n", error, num_nod_vars);
    error = ex_put_var_names (exoid, "n", num_nod_vars, nodal_var_names);
    fprintf (log_file,"after ex_put_var_names, error = %d\n", error);
    
    
    for(i = 0; i < num_time_steps; i++){
      int whole_time_step = i + 1;
      error = ex_put_time (exoid, whole_time_step, time_values+i);
      fprintf (log_file,"after ex_put_time, error = %d\n", error);
      
      
      for (k=1; k<=num_nod_vars; k++)
	{
	  error = ex_put_nodal_var (exoid, whole_time_step, k, num_nodes,
				    nod_var_values[i][k-1]);
	  fprintf (log_file,"after ex_put_nodal_var, error = %d\n", error);
	  if (error) {
	    ex_close (exoid);
	    exit(-1);
	  }
	}
    }
  }



  error = ex_close (exoid2);
  fprintf (log_file,"after ex_close (2), error = %d\n", error);
  return 0;
}

int free_memory(void)
{


  free (xcoor);
  free (ycoor);
  if (num_dim >= 3)
    free (zcoor);

  free (elem_map);

  for(i = 0; i < num_elem_blk; i++)free(connect[i]);


  /* write element block properties */
  if(num_eb_props > 1)
    {
      for (i=0; i<num_eb_props; i++)free(eb_prop_names[i]);
      free (ebids);
    }

  /* write individual node sets */

  for (i=0; i < num_node_sets; i++)
    {
      free(ns_node_list[i]);
    }
  free(ns_node_list);

  free (nsids);
  free(num_df_ns_in_set);
  



  if((num_ns_props > 0) && (num_node_sets > 0))
    {
      for (i=0; i<num_ns_props; i++)
	free(ns_prop_names[i]);
      free(ns_prop_values);
    }




  /* write individual side sets */
  if(num_side_sets > 0)
    {
      for (i=0; i<num_side_sets; i++)
	{
	  if (num_df_ss_in_set[i] > 0) free (ss_dist_fact[i]);
/* 	  free(node_ctr_list[i]); */
	  free(elem_list[i]);
	  free(side_list[i]);
	}
/*       free(ss_dist_fact[i]); */
/*       free(node_ctr_list); */
      free(elem_list);
      free(side_list);
	  
      for (i=0; i<num_ss_props; i++)
	free(ss_prop_names[i]);
      free (ssids);

    }

  if(num_info > 0)
    {

      for (i=0; i<num_info; i++)
	{
	  free(info[i]);
	}
    }

  free_nodal_vars();

  return 0;
}

void free_nodal_vars(void){
  for (i=0; i<num_nod_vars; i++)
    {
      free(nodal_var_names[i]);
    }
  
  if(num_nod_vars > 0){
    free(nodal_var_names);
  }  
  
  if(num_time_steps > 0){
    free(time_values);
  }

  for(i = 0; i < num_time_steps; i++){
    for(j = 0;j < num_nod_vars;j++){
      free(nod_var_values[i][j]);
    }
    free(nod_var_values[i]);
  }

  if(num_time_steps > 0 && num_nod_vars > 0){
    free(nod_var_values);
  }  

  num_nod_vars = 0;
  num_time_steps = 0;

}


void setmem_x_y_z()
{
  xcoor = (double *) calloc(num_nodes, sizeof(double));
  ycoor = (double *) calloc(num_nodes, sizeof(double));
  zcoor = (double *) calloc(num_nodes, sizeof(double));
  return;
}

void copy_resize_x_y_z(int new_nodes)
{
  double * tx,*ty,*tz;
  tx = (double *) calloc(num_nodes+new_nodes, sizeof(double));
  ty = (double *) calloc(num_nodes+new_nodes, sizeof(double));
  tz = (double *) calloc(num_nodes+new_nodes, sizeof(double));

  memset(tx,0,(num_nodes+new_nodes)*sizeof(double));
  memset(ty,0,(num_nodes+new_nodes)*sizeof(double));
  memset(tz,0,(num_nodes+new_nodes)*sizeof(double));

  memcpy(tx,xcoor,num_nodes*sizeof(double));
  memcpy(ty,ycoor,num_nodes*sizeof(double));
  memcpy(tz,zcoor,num_nodes*sizeof(double));

  free(xcoor);
  free(ycoor);
  free(zcoor);
  
  xcoor = tx;
  ycoor = ty;
  zcoor = tz;
  num_nodes += new_nodes;
  return;
}


void setmem_el_blk()
{
  ebids                = (int *) calloc(num_elem_blk, sizeof(int));
  num_elem_in_block    = (int *) calloc(num_elem_blk, sizeof(int));
  num_nodes_per_elem   = (int *) calloc(num_elem_blk, sizeof(int));
  num_attr             = (int *) calloc(num_elem_blk, sizeof(int));
}


void re_setmem_el_blk(int new_size)
{
  int * tebids,*tnum_elem_in_block,*tnum_nodes_per_elem,*tnum_attr;
  int ncpy;
  
  tebids                = (int *) calloc(new_size, sizeof(int));
  tnum_elem_in_block    = (int *) calloc(new_size, sizeof(int));
  tnum_nodes_per_elem   = (int *) calloc(new_size, sizeof(int));
  tnum_attr             = (int *) calloc(new_size, sizeof(int));
  
  memset(tebids,0,new_size*sizeof(int));
  memset(tnum_elem_in_block,0,new_size*sizeof(int));
  memset(tnum_nodes_per_elem,0,new_size*sizeof(int));
  memset(tnum_attr,0,new_size*sizeof(int));
  
  ncpy = num_elem_blk;
  if(new_size < num_elem_blk)ncpy = new_size;
  
  memcpy(tebids,ebids,ncpy*sizeof(int));
  memcpy(tnum_elem_in_block,num_elem_in_block,ncpy*sizeof(int));
  memcpy(tnum_nodes_per_elem,num_nodes_per_elem,ncpy*sizeof(int));
  memcpy(tnum_attr,num_attr,ncpy*sizeof(int));

  free(ebids);
  free(num_elem_in_block);
  free(num_nodes_per_elem);
  free(num_attr);

  ebids                = tebids;
  num_elem_in_block    = tnum_elem_in_block;
  num_nodes_per_elem   = tnum_nodes_per_elem;
  num_attr             = tnum_attr;
  
/*   num_elem_blk = new_size; */
  
}

void setmem_el_connect()
{
  connect              = (int **) calloc(num_elem_blk, sizeof(int*));
  for (i=0; i<num_elem_blk; i++){
    connect[i] = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), sizeof(int));
  }
}

void re_setmem_el_connect(int new_size,int * old_num_elem)
{
  int ** tconnect = (int **) calloc(new_size, sizeof(int*));
  for (i=0; i<new_size; i++){
    tconnect[i] = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]), sizeof(int));
  }

  for(i = 0; i < num_elem_blk; i ++){
    memcpy(tconnect[i],connect[i],num_nodes_per_elem[i] * old_num_elem[i]*sizeof(int));
  }

  for(i = 0; i < num_elem_blk; i ++){
    free(connect[i]);
  }

  free(connect);
  
  connect = tconnect;

}

void setmem_coord_names()
{
  for (i=0; i<num_dim; i++)
    {
      coord_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }
  strcpy(coord_names[0],"XCOOR");
  strcpy(coord_names[1],"YCOOR");
  strcpy(coord_names[2],"ZCOOR");

}

void setmem_elem_map()
{
  int i ;
  elem_map= (int *) calloc(num_elem, sizeof(int));
  for( i = 0; i < num_elem; i = i + 1)
    {
      elem_map[i] = i;
    }
}


void re_setmem_elem_map()
{
  int i ;
  free(elem_map);
  elem_map= (int *) calloc(num_elem, sizeof(int));
  for( i = 0; i < num_elem; i = i + 1)
    {
      elem_map[i] = i;
    }
}


void re_setmem_ns()
{

  free(nsids);
  free(num_nodes_in_set);
  free(num_df_ns_in_set);
  free(ns_node_list);
  free(ns_prop_values);

  nsids = (int *) calloc(num_node_sets, sizeof(int));
  num_nodes_in_set   = (int *) calloc(num_node_sets, sizeof(int));
  num_df_ns_in_set      = (int *) calloc(num_node_sets, sizeof(int));
  ns_node_list       = (int **) calloc(num_node_sets,sizeof(int*));
  ns_prop_values = (int *) calloc (num_node_sets, sizeof(int));

}
void setmem_ns()
{
  nsids = (int *) calloc(num_node_sets, sizeof(int));
  num_nodes_in_set   = (int *) calloc(num_node_sets, sizeof(int));
  num_df_ns_in_set      = (int *) calloc(num_node_sets, sizeof(int));
  ns_node_list       = (int **) calloc(num_node_sets,sizeof(int*));
}


void setmem_ss()
{
  ssids = (int *) calloc(num_side_sets, sizeof(int));
  num_sides_in_set = (int *) calloc(num_side_sets, sizeof(int));
  num_df_ss_in_set = (int *) calloc(num_side_sets, sizeof(int));
  num_elem_in_set = (int *) calloc(num_side_sets, sizeof(int));
  elem_list    = (int **) calloc(num_side_sets,sizeof(int*));
  side_list    = (int **) calloc(num_side_sets,sizeof(int*));
}

int * reset_int_star(int old_size, int added_size,int * orig){
  int * t;
  t = (int *) calloc(old_size+added_size, sizeof(int));
  memcpy(t,orig,old_size*sizeof(int));
  return t;
}

void re_setmem_ss(int added_number)
{
  int ict;
  int new_size = num_side_sets+added_number;
  int *tssids = reset_int_star(num_side_sets,added_number,ssids);
  free(ssids);
  ssids = tssids;
  int *tnum_sides_in_set = reset_int_star(num_side_sets,added_number,num_sides_in_set);
  free(num_sides_in_set);
  num_sides_in_set = tnum_sides_in_set;
  int *tnum_df_ss_in_set = reset_int_star(num_side_sets,added_number,num_df_ss_in_set);
  free(num_df_ss_in_set);
  num_df_ss_in_set = tnum_df_ss_in_set;
  int *tnum_elem_in_set = reset_int_star(num_side_sets,added_number,num_elem_in_set);
  free(num_elem_in_set);
  num_elem_in_set = tnum_elem_in_set;
  int **telem_list    = (int **) calloc(new_size,sizeof(int*));
  int **tside_list    = (int **) calloc(new_size,sizeof(int*));
  double **tss_dist_fact = (double **) calloc(new_size,sizeof(double*));

  for(ict = 0; ict < new_size; ict ++){
    telem_list[ict] = NULL;
    tside_list[ict] = NULL;
    tss_dist_fact[ict] = NULL;
  }

  for(ict = 0; ict < num_side_sets; ict ++){
    telem_list[ict] = elem_list[ict];
    tside_list[ict] = side_list[ict];
    tss_dist_fact[ict] = ss_dist_fact[ict];
  }
  free(elem_list);
  free(side_list);
  free(ss_dist_fact);
  elem_list = telem_list;
  side_list = tside_list;
  ss_dist_fact = tss_dist_fact;

 

}


