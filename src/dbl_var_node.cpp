/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "dbl_var_node.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

DblVarNode::DblVarNode() : VarNode()
{
  type = VAR_DBL;
}

/* ---------------------------------------------------------------------- */
// double DblVarNode::go(int index_in)
// {
//   return data_p[index_in];
// }

/* ---------------------------------------------------------------------- */

void DblVarNode::set_data_pointer(void *in)
{
  data_p = (double *)(in);
}
