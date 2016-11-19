/*
// Copyright (c) 2016 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <m_pd.h>

static t_class *vbap_class;
static t_symbol* s_vbap_2d;
static t_symbol* s_vbap_3d;

typedef struct _vbap_pt
{
    t_float x;
    t_float y;
    t_float z;
} t_vbap_pt;

typedef struct _vbap
{
    t_object    v_obj;
    t_vbap_pt*  v_pts;
    size_t      v_npts;
    
    t_vbap_pt*  v_pts_selected;
    size_t      v_npts_selected;
} t_vbap;

static void *vbap_new(t_symbol *s, int argc, t_atom *argv)
{
    t_vbap_pt* x = (t_vbap_pt*)pd_new(vbap_class);
    if(x)
    {
        
    }
    return x;
}

static void vbap_bang(t_object *x)
{
    
}

extern void vbap_setup(void)
{
    vbap_class = class_new(gensym("vbap"), (t_newmethod)vbap_new, (t_method)NULL, sizeof(t_vbap), CLASS_DEFAULT, A_GIMME, 0);
    class_addbang(vbap_class, (t_method)vbap_bang);
    
    s_vbap_2d = gensym("2d");
    s_vbap_3d = gensym("3d");
}
