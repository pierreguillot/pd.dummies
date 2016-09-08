/*
 // Copyright (c) 2016 Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include <m_pd.h>

static t_class *kop_class;

typedef struct _kop
{
    t_object    m_obj;
    t_float     m_f;
    t_float     m_id;
} t_kop;

static void *kop_new(float f)
{
    t_kop *x = (t_kop *)pd_new(kop_class);
    if(x)
    {
        x->m_id = f;
        signalinlet_new((t_object *)x, 0);        
        outlet_new((t_object *)x, &s_signal);
    }
    return (x);
}

void kop_dsp(t_kop *x, t_signal **sp)
{
    post("kop %i : [%lx %lx] [%lx]", (int)x->m_id, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec);
}

extern void kop_tilde_setup(void)
{
    t_class* c = class_new(gensym("kop~"), (t_newmethod)kop_new, NULL, sizeof(t_kop), CLASS_DEFAULT, A_FLOAT, 0);
    if(c)
    {
        class_addmethod(c, (t_method)kop_dsp, gensym("dsp"), A_CANT);
        CLASS_MAINSIGNALIN(c, t_kop, m_f);
    }
    kop_class = c;
}
