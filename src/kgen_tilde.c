/*
 // Copyright (c) 2016 Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include <m_pd.h>

static t_class *kgen_class;

typedef struct _kgen
{
    t_object    m_obj;
    t_float     m_f;
    t_float     m_id;
} t_kgen;

static void *kgen_new(float f)
{
    t_kgen *x = (t_kgen *)pd_new(kgen_class);
    if(x)
    {
        x->m_id = f;
        outlet_new((t_object *)x, &s_signal);
    }
    return (x);
}

void kgen_dsp(t_kgen *x, t_signal **sp)
{
    post("kgen %i : [%lx]", (int)x->m_id, sp[0]->s_vec);
}

extern void kgen_tilde_setup(void)
{
    t_class* c = class_new(gensym("kgen~"), (t_newmethod)kgen_new, NULL, sizeof(t_kgen), CLASS_NOINLET, A_FLOAT, 0);
    if(c)
    {
        class_addmethod(c, (t_method)kgen_dsp, gensym("dsp"), A_CANT);
    }
    kgen_class = c;
}
