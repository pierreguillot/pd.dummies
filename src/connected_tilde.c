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
    }
    return x;
}

static void connected_tilde_dsp(t_connected *x, t_signal **sp)
{
    int i;
    t_linetraverser t;
    t_outconnect*   oc;
    char*           st;
    
    st = getbytes(sizeof(char) * (size_t)obj_ninlets((t_object *)x));
    if(st)
    {
        linetraverser_start(&t, x->c_cnv);
        while((oc = linetraverser_next(&t)))
        {
            if(t.tr_ob2 == (t_object *)x && obj_issignaloutlet(t.tr_ob, t.tr_outno))
            {
                st[t.tr_inno] = (char)1;
            }
        }
        startpost("connections %ld : ", (unsigned int)x);
        for(i = 0; i < obj_ninlets((t_object *)x); ++i)
        {
            postfloat((float)st[i]);
            poststring(" ");
        }
        endpost();
        freebytes(st, sizeof(char) * (size_t)obj_ninlets((t_object *)x));
    }
    else
    {
        pd_error(x, "can't allocate memory for inputs states.");
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
