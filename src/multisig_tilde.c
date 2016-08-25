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
    t_sample*   m_temp;
    size_t      m_size;
} t_multi;

static void *multi_new(t_symbol *s, int argc, t_atom *argv)
{
    t_multi *x = (t_multi *)pd_new(multi_class);
    if(x)
    {
        signalinlet_new((t_object *)x, 0);
        signalinlet_new((t_object *)x, 0);
        signalinlet_new((t_object *)x, 0);
        
        outlet_new((t_object *)x, &s_signal);
        outlet_new((t_object *)x, &s_signal);
        outlet_new((t_object *)x, &s_signal);
        outlet_new((t_object *)x, &s_signal);
        
        x->m_size = 0;
        x->m_temp = NULL;
    }
    
    return (x);
}

void multi_free(t_multi *x)
{
    if(x->m_size && x->m_temp)
    {
        freebytes(x->m_temp, x->m_size);
        x->m_temp = NULL;
        x->m_size = 0;
    }
}

t_int *multi_perform(t_int *w)
{
    int i;
    t_multi   *x   = (t_multi *)(w[1]);
    t_sample  *in1 = (t_sample *)(w[2]);
    t_sample  *in2 = (t_sample *)(w[3]);
    t_sample  *in3 = (t_sample *)(w[4]);
    t_sample  *in4 = (t_sample *)(w[5]);
    t_sample  *out1 = (t_sample *)(w[6]);
    t_sample  *out2 = (t_sample *)(w[7]);
    t_sample  *out3 = (t_sample *)(w[8]);
    t_sample  *out4 = (t_sample *)(w[9]);
    int          n = (int)(w[10]);
    t_sample  *temp = x->m_temp;
    
    for(i = 0; i < n; ++i) {
        *temp++ = *in1++;
    }
    for(i = 0; i < n; ++i) {
        *temp++ = *in2++;
    }
    for(i = 0; i < n; ++i) {
        *temp++ = *in3++;
    }
    for(i = 0; i < n; ++i) {
        *temp++ = *in4++;
    }
    
    temp = x->m_temp;
    for(i = 0; i < n; ++i) {
        *out1++ = *temp++;
    }
    for(i = 0; i < n; ++i) {
        *out2++ = *temp++;
    }
    for(i = 0; i < n; ++i) {
        *out3++ = *temp++;
    }
    for(i = 0; i < n; ++i) {
        *out4++ = *temp++;
    }
    
    return (w+11);
}

void multi_dsp(t_multi *x, t_signal **sp)
{
    post("ins : %lx %lx %lx %lx", sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec);
    post("outs : %lx %lx %lx %lx", sp[4]->s_vec, sp[5]->s_vec, sp[6]->s_vec, sp[7]->s_vec);
    
    if(x->m_size && x->m_temp)
    {
        freebytes(x->m_temp, x->m_size);
        x->m_temp = NULL;
        x->m_size = 0;
    }
    x->m_size = sizeof(t_sample) * 4 * sp[0]->s_n;
    if(x->m_size)
    {
        x->m_temp = (t_sample *)getbytes(x->m_size);
        if(x->m_temp)
        {
            dsp_add(multi_perform, 10, x,
                    sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec,
                    sp[4]->s_vec, sp[5]->s_vec, sp[6]->s_vec, sp[7]->s_vec, sp[0]->s_n);
        }
        else
        {
            pd_error(x, "can't allocate temporary vectors.");
        }
    }
    
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
