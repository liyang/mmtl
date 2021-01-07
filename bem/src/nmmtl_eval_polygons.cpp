
/*
  
  FACILITY:  NMMTL
  
  MODULE DESCRIPTION:
  
  Contains nmmtl_evaluate_polygons function.
  
  AUTHORS:
  
  Kevin J. Buchs
  
  CREATION DATE:  Fri Nov 22 10:40:38 1991
  
  COPYRIGHT:   Copyright (C) 1991 by Mayo Foundation. All rights reserved.
  
  
  */


/*
 *******************************************************************
 **  INCLUDE FILES
 *******************************************************************
 */

#include "nmmtl.h"


/*
 *******************************************************************
 **  STRUCTURE AND TYPE DEFINITIONS
 *******************************************************************
 */

/* structure to hold info about the polygon points */

typedef struct pgnpts
{
  double x[2],y[2];         /* the coordinates of start [0] and end [1] */
                            /* of the vector */
  double dx,dy;             /* the displacements of the vector in x and */
                            /* y */
  float theta2[2];          /* the edge angles */
  int valid;                /* indicates if this entry on list is used */
  struct pgnpts *next;
} PGNPTS, *PGNPTS_P;




/*
 *******************************************************************
 **  MACRO DEFINITIONS
 *******************************************************************
 */
/*
 *******************************************************************
 **  PREPROCESSOR CONSTANTS
 *******************************************************************
 */
/*
 *******************************************************************
 **  GLOBALS
 *******************************************************************
 */

#ifdef NMMTL_DUMP_DIAG
extern FILE *dump_file;  /* a file for diagnostics */
#endif

/*
 *******************************************************************
 **  FUNCTION DECLARATIONS
 *******************************************************************
 */
/*
 *******************************************************************
 **  FUNCTION DEFINITIONS
 *******************************************************************
 */


/*
  
  FUNCTION NAME:  nmmtl_evaluate_polygons
  
  FUNCTIONAL DESCRIPTION:
  
  Takes a polygon contour and processes it into a series of line
  segments, with appropriate angle information.
  
  
  FORMAL PARAMETERS:
  
  int cntr_seg
  number of segments to break a contour into
  
  float half_minimum_dimension
  half of the smallest geometric dimension - used to determine if
  segments are broken small enough.
  
  int conductor_counter
  keep count of the conductors
  
  CONTOURS_P contour:
  Conductor or groundwire contour.
  
  LINE_SEGMENTS_P *segments:
  output list of line segments generated by fracturing the 
  polygon found.
  
  EXTENT_DATA_P extent_data:
  structure of extents of conductor region

  RETURN VALUE:
  
  SUCCESS or FAIL
  
  DESIGN:
  
  
  As we look at a polygon, we must traverse it generating vectors based
  on each pair of verticies.  We need to determine if each vertex is
  convex (points out) or concave (points in).  This can be done by 
  summing the angles turned as the perimeter of the polygon is traversed
  where a left hand turn is positive and a right hand is negative.
  The sum of all the turns should be +360 degrees or -360 degrees, 
  depending upon which direction the polygon is traversed.  When this
  sum has been computed the sign will allow us to determine if each
  vertex is convex (sign matches sign of sum) or concave (sign opposite
  the sign of the sum).
  
  We do not need to further consider concave verticies.  For the convex
  verticies, we need to compute the angle swung by an arc on the outside
  of the polygon.  This will be the angle turned plus 180 degrees.
  
  Theta2[0] is the angle formed by the edge containing the leading end
  of the vector, and theta2[1] is the other edge.  If one of the ends of
  the vector does not form an convex edge, then the edge angle is zero.
  
  CALLING SEQUENCE:
  
  status = nmmtl_evaluate_polygons(cntr_seg,half_minimum_dimension,
  conductor_counter,contours,&segments,extent_data)
  
  
  */


int nmmtl_evaluate_polygons(int cntr_seg,
#ifndef NO_HALF_MIN_CHECKING
					float half_minimum_dimension,
#endif
					int conductor_counter,
					CONTOURS_P contour,
					LINE_SEGMENTS_P *segments,
					EXTENT_DATA_P extent_data)
{
  static PGNPTS_P head = NULL,last,current;
  POLYPOINTS_P point, last_point;
  int i;
  float sum_of_angles;
  int sign_of_polygon, sign_of_angle;
  LINE_SEGMENTS_P new_segment, last_segment, leading_segment = NULL;
  
  
  if(contour->points == NULL || contour->points->next == NULL)
    return(FAIL);
  
  /* preallocate a linked list that can remain and be extended if need be */
  
  if(head == NULL)
  {
    head = (PGNPTS_P)malloc(sizeof(PGNPTS));
    head->valid = 0;
    current = head;
    for(i=1;i < 10; i++)
    {
      current->next = (PGNPTS_P)malloc(sizeof(PGNPTS));
      current = current->next;
      current->valid = 0;
    }
    current->next = NULL;
    current = head;
  }
  
  /* if list was preallocated, mark it all invalid */
  else
  {
    for(current=head; current != NULL; current = current->next)
    {
      current->valid = 0;
    }
  }
  
  /* find the end of the segment list */
  last_segment = *segments;
  if(last_segment != NULL)
  {
    while(last_segment->next != NULL)
      last_segment = last_segment->next;
  }
  
  
  
  
  last_point = contour->points;

  /* set up the new extent data */

  /* left conductor extent */
  if(extent_data->left_cond_extent > last_point->x)
    extent_data->left_cond_extent = last_point->x;

  /* right conductor extent */
  if(extent_data->right_cond_extent < last_point->x)
    extent_data->right_cond_extent = last_point->x;

  /* minimum conductor height above ground plane */
  if(extent_data->min_cond_height > last_point->y)
    extent_data->min_cond_height = last_point->y;

  point = contour->points->next;
  current = head;
  current->valid = 0;
  
  while(point != NULL)
  {
    current->valid = 1;
    current->x[0] = last_point->x;
    current->x[1] = point->x;
    current->dx = point->x - last_point->x;
    current->y[0] = last_point->y;
    current->y[1] = point->y;
    current->dy = point->y - last_point->y;

    /* set up the new extent data */

    /* left conductor extent */
    if(extent_data->left_cond_extent > point->x)
      extent_data->left_cond_extent = point->x;

    /* right conductor extent */
    if(extent_data->right_cond_extent < point->x)
      extent_data->right_cond_extent = point->x;

    /* minimum conductor height above ground plane */
    if(extent_data->min_cond_height > point->y)
      extent_data->min_cond_height = point->y;


    
    /* do we need a longer list ? */
    if(current->next == NULL)
    {
      last = current;
      
      /* add ten more entries */
      for(i=1;i < 10; i++)
      {
	      last->next = (PGNPTS_P)malloc(sizeof(PGNPTS));
	      last = last->next;
	      last->valid = 0;
      }
      last->next = NULL;
    }
    
    current = current->next;
    current->valid = 0;
    last_point = point;
    point = point->next;
  }
  
  /* now compute the angles */
  sum_of_angles = 0;
  last = head;
  current = head->next;
  do
  {
    last->theta2[1] = nmmtl_angle_of_intersection(last->dx,last->dy,
						  current->dx,current->dy);
    current->theta2[0] = last->theta2[1];
    sum_of_angles += last->theta2[1];
    last = current;
    current = current->next;
    
  } while(current->valid != 0);
  
  /* and then add in the last turn around to the first vector */
  
  current = head;
  last->theta2[1] = nmmtl_angle_of_intersection(last->dx,last->dy,
						current->dx,current->dy);
  current->theta2[0] = last->theta2[1];
  sum_of_angles += last->theta2[1];
  
#ifdef DIAG_POLY
#ifdef NMMTL_DUMP_DIAG
  fprintf(dump_file,"polygon vector:\n");
  for(current=head; current->valid != 0; current = current->next)
  {
    fprintf(dump_file,"  (%f,%f)   angle %f degrees\n",
	    current->dx,current->dy,
	    current->theta2[1] * RADIANS_TO_DEGREES );
  }
  fprintf(dump_file,"sum of angles: %f\n",
	  sum_of_angles * RADIANS_TO_DEGREES);
#endif
#endif
  
  /* have now determined the sum of the angles turned when going */
  /* around the polygon */
  
  sign_of_polygon = sum_of_angles < 0 ? -1 : 1;
  
  /* throw away concave verticies, add 180 degrees to convex ones to */
  /* get the true exterior arc angle */
  
  for(current=head; current->valid != 0; current = current->next)
  {
    /* adjust theta2[0] */
    sign_of_angle = current->theta2[0] < 0 ? -1 : 1;
    if(sign_of_angle == sign_of_polygon)
      current->theta2[0] = PI + fabs(current->theta2[0]);
    else
      current->theta2[0] = 0;
    
    /* adjust theta2[1] */
    sign_of_angle = current->theta2[1] < 0 ? -1 : 1;
    if(sign_of_angle == sign_of_polygon)
      current->theta2[1] = PI + fabs(current->theta2[1]);
    else
      current->theta2[1] = 0;
  }	
  
#ifdef DIAG_POLY
#ifdef NMMTL_DUMP_DIAG
  for(current=head; current->valid != 0; current = current->next)
  {
    fprintf(dump_file,"  (%f,%f)   theta2 %f degrees\n",
	    current->dx,current->dy,
	    current->theta2[1] * RADIANS_TO_DEGREES );
  }
#endif
#endif
  
  for(current=head; current->valid != 0; current = current->next)
  {
    /* for each edge of the polygon, set up a new line segment */
    /* with the appropriate parameters */
    
    new_segment = (LINE_SEGMENTS_P)malloc(sizeof(LINE_SEGMENTS));
    new_segment->startx = current->x[0];
    new_segment->endx = current->x[1];
    new_segment->starty = current->y[0];
    new_segment->endy = current->y[1];
    new_segment->theta2[0] = current->theta2[0];
    new_segment->theta2[1] = current->theta2[1];
    
    if(current->theta2[0] != 0)
    {
      new_segment->nu[0] = PI/current->theta2[0];
      new_segment->free_space_nu[0] = PI/current->theta2[0];
    }
    else
    {
      new_segment->nu[0] = 0;
      new_segment->free_space_nu[0] = 0;
    }
    if(current->theta2[1] != 0)
    {
      new_segment->nu[1] = PI/current->theta2[1];
      new_segment->free_space_nu[1] = PI/current->theta2[1];
    }
    else
    {
      new_segment->nu[1] = 0;
      new_segment->free_space_nu[1] = 0;
    }
    
    new_segment->conductor = conductor_counter;
    new_segment->epsilon[0] = 0.0;
    new_segment->epsilon[1] = 0.0;
    if(leading_segment == NULL) leading_segment = new_segment;
    
    /* which way do you turn to get to the interior of polygon - same */
    /* direction that the polygon turns */
    new_segment->interior = sign_of_polygon;
    
    new_segment->length = 
      sqrt(current->dx * current->dx + current->dy * current->dy);
    
    /* number of elements this should be broken into - take the */
    /* ratio of the length of the side to the total perimeter of the */
    /* polygon and multiply by the number of segments to break each */
    /* contour into.  contour->x0 is the perimeter of the polygon */
    /* under consideration.  Round up to the next higher integer */
    
    new_segment->divisions = (int)(.99 + cntr_seg * new_segment->length /
      contour->x0);
    
    /* force two elements per side, minimum */
    if(new_segment->divisions < 2) new_segment->divisions = 2;

#ifndef NO_HALF_MIN_CHECKING
    
    /* is this element size too large?  Check to see if it is greater */
    /* than half of the smallest conductor edge. */
    
    if(new_segment->length/new_segment->divisions > half_minimum_dimension)
    {
      new_segment->divisions = .99 + new_segment->length / 
	half_minimum_dimension;
    }
    
#endif    
    
    /* associate edge pairs */
    /* for first segment, edge pair association will need to */
    /* wait till the list is done */
    if(leading_segment != new_segment)
    {
      last_segment->edge_pair[1] = new_segment;
      new_segment->edge_pair[0] = last_segment;
    }
    
    /* hook up to end of the list */
    if(last_segment == NULL)
    {
      *segments = new_segment;
      last_segment = new_segment;
    }
    else
    {
      last_segment->next = new_segment;
      last_segment = new_segment;
    }
  }
  
  /* final segment is complete : */
  
  /* terminate list */
  new_segment->next = NULL;
  
  /* now make last edge pair assocation */
  new_segment->edge_pair[1] = leading_segment;
  leading_segment->edge_pair[0] = new_segment;
  
  return(SUCCESS);
}




