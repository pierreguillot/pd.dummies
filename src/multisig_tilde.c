/*
// Copyright (c) 2016 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <m_pd.h>

static t_class *multi_class;

typedef struct _multi
{
    t_object    m_obj;
    t_float     m_f;
} t_multi;

static void *multi_new(t_symbol *s, int argc, t_atom *argv)
{
    t_object *x = (t_object *)pd_new(multi_class);
    if(x)
    {
        signalinlet_new(x, 0);
        signalinlet_new(x, 0);
        signalinlet_new(x, 0);
        
        outlet_new(x, &s_signal);
        outlet_new(x, &s_signal);
        outlet_new(x, &s_signal);
        outlet_new(x, &s_signal);
    }
    
    return (x);
}

t_int *multi_perform(t_int *w)
{
    int i;
    t_sample  *in1 = (t_sample *)(w[2]);
    t_sample  *in2 = (t_sample *)(w[3]);
    t_sample  *in3 = (t_sample *)(w[4]);
    t_sample  *in4 = (t_sample *)(w[5]);
    t_sample  *out1 = (t_sample *)(w[6]);
    t_sample  *out2 = (t_sample *)(w[7]);
    t_sample  *out3 = (t_sample *)(w[8]);
    t_sample  *out4 = (t_sample *)(w[9]);
    int          n = (int)(w[10]);
    
    for(i = 0; i < n; ++i)
    {
        *out1++ = *in1++;
        *out2++ = *in2++;
        *out3++ = *in3++;
        *out4++ = *in4++;
    }
    return (w+11);
}

void multi_dsp(t_object *x, t_signal **sp)
{
    post("ins : %lx %lx %lx %lx", sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec);
    post("outs : %lx %lx %lx %lx", sp[4]->s_vec, sp[5]->s_vec, sp[6]->s_vec, sp[7]->s_vec);
    dsp_add(multi_perform, 10, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec,
            sp[4]->s_vec, sp[5]->s_vec, sp[6]->s_vec, sp[7]->s_vec, sp[0]->s_n);
    
}


extern void multisig_tilde_setup(void)
{
    t_class* c = class_new(gensym("multisig~"), (t_newmethod)multi_new, NULL, sizeof(t_multi), CLASS_DEFAULT, 0);
    if(c)
    {
        class_addmethod(c, (t_method)multi_dsp, gensym("dsp"), A_CANT);
        CLASS_MAINSIGNALIN(c, t_multi, m_f);
    }
    multi_class = c;
}
