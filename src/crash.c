/*
// Copyright (c) 2016 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <m_pd.h>

static t_class *crash_class;

static void *crash_new(t_symbol *s, int argc, t_atom *argv)
{
    return pd_new(crash_class);
}

static void crash_bang(t_object *x)
{
    sys_ouch();
}

extern void crash_setup(void)
{
    crash_class = class_new(gensym("crash"), (t_newmethod)crash_new, (t_method)NULL, sizeof(t_object), CLASS_DEFAULT, 0);
    class_addbang(crash_class, (t_method)crash_bang);
}
