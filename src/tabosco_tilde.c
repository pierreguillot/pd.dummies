/*
 // Copyright (c) 2016 Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include <m_pd.h>
#include <g_canvas.h>
#include <m_imp.h>

typedef struct _tabosco
{
    t_object    t_obj;
    t_symbol*   t_name;
    t_word*     t_buffer;
    int         t_size;
    t_float     t_dummy;
    t_sample    t_index;
} t_tabosco;

static t_class *tabosco_tilde_class;

static void *tabosco_tilde_new(t_symbol* s)
{
    t_tabosco* x = (t_tabosco *)pd_new(tabosco_tilde_class);
    if(x)
    {
        x->t_name   = s;
        x->t_index  = 0.f;
        outlet_new((t_object *)x, &s_signal);
    }
    return x;
}

t_int *tabosco_perform(t_int *w)
{
    t_tabosco* x    = (t_tabosco *)(w[1]);
    t_sample*  in   = (t_sample *)(w[2]);
    t_sample*  out  = (t_sample *)(w[3]);
    size_t n        = (size_t)(w[4]);
    t_sample sr     = (t_sample)(w[5]);
    int size        = x->t_size;
    t_word* buffer  = x->t_buffer;
    t_sample index  = x->t_index;
    t_sample freq;
    while(n--)
    {
        freq = *in++ / sr;
        index += freq;
        if(freq < 0.f && index < 0.f) {
            index += 1.f;
        }
        else if(freq > 0.f && index >= 1.f) {
            index -= 1.f;
        }
        
        *out++ = buffer[(int)(index * (t_sample)size)].w_float;
    }
    x->t_index = index;
    
    return (w+6);
}


static void tabosco_tilde_dsp(t_tabosco *x, t_signal **sp)
{
    t_garray *a = NULL;
    if(!x->t_name)
    {
        pd_error(x, "tabosco~: ain't got no name.");
        return;
    }
    if (!(a = (t_garray *)pd_findbyclass(x->t_name, garray_class)))
    {
        pd_error(x, "tabosco~: %s no such array.", x->t_name->s_name);
        return;
    }
    else if (!garray_getfloatwords(a, &x->t_size, &x->t_buffer))
    {
        pd_error(x, "tabosco~: %s array is empty.", x->t_name->s_name);
        return;
    }
    else
    {
        garray_usedindsp(a);
    }
    dsp_add(tabosco_perform, 5, (t_int)x, (t_int)sp[0]->s_vec, (t_int)sp[1]->s_vec, (t_int)sp[0]->s_n, (t_int)sp[0]->s_sr);
}

void tabosco_tilde_setup(void)
{
    t_class* c = class_new(gensym("tabosco~"), (t_newmethod)tabosco_tilde_new, (t_method)NULL, sizeof(t_tabosco), CLASS_DEFAULT, A_SYMBOL, 0);
    if(c)
    {
        class_addmethod(c, (t_method)tabosco_tilde_dsp, gensym("dsp"), A_CANT);
        CLASS_MAINSIGNALIN(c, t_tabosco, t_dummy);
    }
    tabosco_tilde_class = c;
}
