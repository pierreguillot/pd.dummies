/*
// Copyright (c) 2016 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <m_pd.h>
#include <g_canvas.h>
#include <m_imp.h>

typedef struct _connected
{
    t_object    c_obj;
    t_canvas*   c_cnv;
    t_outlet*   c_out;
} t_connected;

static t_class *connected_tilde_class;

static void *connected_tilde_new(float f)
{
    int i;
    t_connected* x = (t_connected *)pd_new(connected_tilde_class);
    if(x)
    {
        x->c_cnv = canvas_getcurrent();
        for(i = 0; i < (int)f; ++i)
        {
            signalinlet_new((t_object *)x, 0);
        }
        x->c_out = outlet_new((t_object *)x, &s_list);
    }
    return x;
}

static void connected_tilde_dsp(t_connected *x, t_signal **sp)
{
    int i, ninlets = obj_ninlets((t_object *)x);
    t_linetraverser t;
    t_outconnect*   oc;
    t_atom*         av = (t_atom *)getbytes(ninlets * sizeof(t_atom));
    if(av)
    {
        for(i = 0; i < ninlets; ++i)
        {
            SETFLOAT(av+i, 0);
        }
        
        linetraverser_start(&t, x->c_cnv);
        while((oc = linetraverser_next(&t)))
        {
            if(t.tr_ob2 == (t_object *)x && obj_issignaloutlet(t.tr_ob, t.tr_outno))
            {
                SETFLOAT(av+t.tr_inno, 1);
            }
        }
        outlet_list(x->c_out, &s_list, ninlets, av);
        freebytes(av, ninlets * sizeof(t_atom));
    }
    else
    {
        pd_error(x, "can't allocate memory to output the connections states.");
    }
}

extern void connected_tilde_setup(void)
{
    t_class* c = class_new(gensym("connected~"), (t_newmethod)connected_tilde_new, (t_method)NULL, sizeof(t_connected), CLASS_NOINLET, A_FLOAT, 0);
    if(c)
    {
        class_addmethod(c, (t_method)connected_tilde_dsp, gensym("dsp"), A_CANT);
    }
    connected_tilde_class = c;
}
