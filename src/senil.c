/*
// Copyright (c) 2016 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <m_pd.h>

static t_class *senil_class;

static void *senil_new(t_symbol *s, int argc, t_atom *argv)
{
    return pd_new(senil_class);
}

static void senil_bang(t_object *x)
{
    sys_ouch();
}

static void senil_signal(t_object *x)
{
    post("senil_signal");
}

static void senil_dsp(t_object *x, t_signal **sp)
{
    post("senil_dsp");
}



extern void senil_setup(void)
{
    senil_class = class_new(gensym("senil"), (t_newmethod)senil_new, (t_method)NULL, sizeof(t_object), CLASS_DEFAULT, 0);
    class_addbang(senil_class, (t_method)senil_bang);
    class_addmethod(senil_class, (t_method)senil_signal, gensym("signal"), A_NULL);
    class_addmethod(senil_class, (t_method)senil_dsp, gensym("dsp"), A_CANT);
}
