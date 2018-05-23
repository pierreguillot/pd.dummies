/*
 // Copyright (c) 2016 Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include <m_pd.h>
#include <g_canvas.h>
#include <m_imp.h>
#include <float.h>
#include <math.h>

typedef struct _maxicatch
{
    t_object    m_obj;
    t_clock*    m_clock;
    t_outlet*   m_outlet;
    t_atom      m_values[2];
    t_float     m_dummy;
} t_maxicatch;

static t_class *maxicatch_tilde_class;

static void maxicatch_tilde_tick(t_maxicatch *x)
{
    outlet_list(x->m_outlet, &s_list, 2, x->m_values);
}

static void *maxicatch_tilde_new()
{
    t_maxicatch* x = (t_maxicatch *)pd_new(maxicatch_tilde_class);
    if(x)
    {
        x->m_clock = clock_new(x, (t_method)maxicatch_tilde_tick);
        x->m_outlet= outlet_new((t_object *)x, &s_list);
    }
    return x;
}

t_int *maxicatch_perform(t_int *w)
{
    t_maxicatch* x    = (t_maxicatch *)(w[1]);
    t_sample*  in   = (t_sample *)(w[2]);
    size_t n        = (size_t)(w[3]);
    int i = 0, index = -1;
    t_sample val = 0.f, lval = 0.f, nval = 0.f, aval = 0.f;

    while(n--)
    {
        nval = *in++;
        aval = fabsf(nval);
        if(aval > lval) {
            val = nval;
            lval = aval;
            index = i;
        }
        i++;
    }
    if(index >= 0)
    {
        SETFLOAT(x->m_values, (float)index);
        SETFLOAT(x->m_values+1, (float)val);
        clock_delay(x->m_clock, 0.f);
    }
    
    return (w+4);
}


static void maxicatch_tilde_dsp(t_maxicatch *x, t_signal **sp)
{
    dsp_add(maxicatch_perform, 3, (t_int)x, (t_int)sp[0]->s_vec, (t_int)sp[0]->s_n);
}

void maxicatch_tilde_setup(void)
{
    t_class* c = class_new(gensym("maxicatch~"), (t_newmethod)maxicatch_tilde_new, (t_method)NULL, sizeof(t_maxicatch), CLASS_DEFAULT, 0);
    if(c)
    {
        class_addmethod(c, (t_method)maxicatch_tilde_dsp, gensym("dsp"), A_CANT);
        CLASS_MAINSIGNALIN(c, t_maxicatch, m_dummy);
    }
    maxicatch_tilde_class = c;
}
