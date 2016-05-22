/*
// Copyright (c) 2016 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <m_pd.h>

static t_class *pak_class;
static t_class* pak_inlet_class;

typedef struct _pak_inlet
{
    t_class*    x_pd;
    t_atom*     x_atoms;
    t_int       x_max;
    t_pd*       x_owner;
} t_pak_inlet;

static void pak_inlet_bang(t_pak_inlet *x)
{
    pd_bang(x->x_owner);
}

static void pak_inlet_float(t_pak_inlet *x, float f)
{
    if(x->x_atoms->a_type == A_FLOAT)
    {
        x->x_atoms->a_w.w_float = f;
        pd_bang(x->x_owner);
    }
    else
    {
        pd_error(x, "pak: wrong type (float)");
    }
}

static void pak_inlet_pointer(t_pak_inlet *x, t_gpointer *gp)
{
    if(x->x_atoms->a_type == A_POINTER)
    {
        gpointer_unset(x->x_atoms->a_w.w_gpointer);
        *(x->x_atoms->a_w.w_gpointer) = *gp;
        if(x->x_atoms->a_w.w_gpointer->gp_stub)
        {
            x->x_atoms->a_w.w_gpointer->gp_stub->gs_refcount++;
        }
        pd_bang(x->x_owner);
    }
    else
    {
        pd_error(x, "pak: wrong type (pointer)");
    }
}

static void pak_inlet_symbol(t_pak_inlet *x, t_symbol* s)
{
    if(x->x_atoms->a_type == A_SYMBOL)
    {
        x->x_atoms->a_w.w_symbol = s;
        pd_bang(x->x_owner);
    }
    else
    {
        pd_error(x, "pak: wrong type (symbol)");
    }
}

static void pak_inlet_list(t_pak_inlet *x, t_symbol* s, int argc, t_atom* argv)
{
    int i;
    for(i = 0; i < x->x_max && i < argc; ++i)
    {
        if(argv[i].a_type == A_FLOAT)
        {
            if(x->x_atoms[i].a_type == A_FLOAT)
            {
                x->x_atoms[i].a_w.w_float = argv[i].a_w.w_float;
            }
            else
            {
                pd_error(x, "pak: wrong type (float)");
            }
        }
        else if(argv[i].a_type == A_POINTER)
        {
            if(x->x_atoms[i].a_type == A_POINTER)
            {
                gpointer_unset(x->x_atoms[i+1].a_w.w_gpointer);
                *(x->x_atoms[i+1].a_w.w_gpointer) = *(argv[i].a_w.w_gpointer);
                if(x->x_atoms[i+1].a_w.w_gpointer->gp_stub)
                {
                    x->x_atoms[i+1].a_w.w_gpointer->gp_stub->gs_refcount++;
                }
            }
            else
            {
                pd_error(x, "pak: wrong type (pointer)");
            }
        }
        if(argv[i].a_type == A_SYMBOL)
        {
            if(x->x_atoms[i].a_type == A_SYMBOL)
            {
                x->x_atoms[i].a_w.w_symbol = argv[i].a_w.w_symbol;
            }
            else
            {
                pd_error(x, "pak: wrong type (symbol)");
            }
        }
    }
    pd_bang(x->x_owner);
}

static void pak_inlet_anything(t_pak_inlet *x, t_symbol* s, int argc, t_atom* argv)
{
    int i;
    if(x->x_atoms[0].a_type == A_SYMBOL)
    {
        x->x_atoms[0].a_w.w_symbol = s;
    }
    else
    {
        pd_error(x, "pak: wrong type (symbol)");
    }
    for(i = 0; i < x->x_max-1 && i < argc; ++i)
    {
        if(argv[i].a_type == A_FLOAT)
        {
            if(x->x_atoms[i+1].a_type == A_FLOAT)
            {
                x->x_atoms[i+1].a_w.w_float = argv[i].a_w.w_float;
            }
            else
            {
                pd_error(x, "pak: wrong type (float)");
            }
        }
        else if(argv[i].a_type == A_POINTER)
        {
            if(x->x_atoms[i+1].a_type == A_POINTER)
            {
                gpointer_unset(x->x_atoms[i+1].a_w.w_gpointer);
                *(x->x_atoms[i+1].a_w.w_gpointer) = *(argv[i].a_w.w_gpointer);
                if(x->x_atoms[i+1].a_w.w_gpointer->gp_stub)
                {
                    x->x_atoms[i+1].a_w.w_gpointer->gp_stub->gs_refcount++;
                }
            }
            else
            {
                pd_error(x, "pak: wrong type (pointer)");
            }
        }
        if(argv[i].a_type == A_SYMBOL)
        {
            if(x->x_atoms[i+1].a_type == A_SYMBOL)
            {
                x->x_atoms[i+1].a_w.w_symbol = argv[i].a_w.w_symbol;
            }
            else
            {
                pd_error(x, "pak: wrong type (symbol)");
            }
        }
    }
    pd_bang(x->x_owner);
}


typedef struct _pak
{
    t_object    x_obj;
    t_int       x_n;            /* number of args */
    t_atom*     x_vec;          /* input values */
    t_int       x_nptr;         /* number of pointers */
    t_gpointer* x_gpointer;     /* the pointers */
    t_atom*     x_outvec;       /* space for output values */
} t_pak;

static void *pak_new(t_symbol *s, int argc, t_atom *argv)
{
    int i;
    int nptr = 0;
    t_pak *x = (t_pak *)pd_new(pak_class);
    t_atom defarg[2], *ap, *vec, *vp;
    t_gpointer *gp;
    if(!argc)
    {
        argv = defarg;
        argc = 2;
        SETFLOAT(&defarg[0], 0);
        SETFLOAT(&defarg[1], 0);
    }

    x->x_n = argc;
    vec = x->x_vec = (t_atom *)getbytes(argc * sizeof(*x->x_vec));
    x->x_outvec = (t_atom *)getbytes(argc * sizeof(*x->x_outvec));

    for(i = argc, ap = argv; i--; ap++)
    {
        if(ap->a_type == A_SYMBOL && *ap->a_w.w_symbol->s_name == 'p')
        {
            nptr++;
        }
    }

    gp = x->x_gpointer = (t_gpointer *)getbytes(nptr * sizeof (*gp));
    x->x_nptr = nptr;

    for (i = 0, vp = x->x_vec, ap = argv; i < argc; i++, ap++, vp++)
    {
        if (ap->a_type == A_FLOAT)
        {
            *vp = *ap;
            if (i) floatinlet_new(&x->x_obj, &vp->a_w.w_float);
        }
        else if (ap->a_type == A_SYMBOL)
        {
            char c = *ap->a_w.w_symbol->s_name;
            if (c == 's')
            {
                SETSYMBOL(vp, &s_symbol);
                if (i) symbolinlet_new(&x->x_obj, &vp->a_w.w_symbol);
            }
            else if (c == 'p')
            {
                vp->a_type = A_POINTER;
                vp->a_w.w_gpointer = gp;
                gpointer_init(gp);
                if (i) pointerinlet_new(&x->x_obj, gp);
                gp++;
            }
            else
            {
                if (c != 'f') pd_error(x, "pak: %s: bad type",
                    ap->a_w.w_symbol->s_name);
                SETFLOAT(vp, 0);
                if (i) floatinlet_new(&x->x_obj, &vp->a_w.w_float);
            }
        }
    }
    outlet_new(&x->x_obj, &s_list);
    return (x);
}

static void pak_bang(t_pak *x)
{
    int i, reentered = 0, size = x->x_n * sizeof (t_atom);
    t_gpointer *gp;
    t_atom *outvec;
    for (i = x->x_nptr, gp = x->x_gpointer; i--; gp++)
        if (!gpointer_check(gp, 1))
    {
        pd_error(x, "pak: stale pointer");
        return;
    }
        /* reentrancy protection.  The first time through use the pre-allocated
        x_outvec; if we're reentered we have to allocate new memory. */
    if (!x->x_outvec)
    {
            /* LATER figure out how to deal with reentrancy and pointers... */
        if (x->x_nptr)
            post("pak_bang: warning: reentry with pointers unprotected");
        outvec = t_getbytes(size);
        reentered = 1;
    }
    else
    {
        outvec = x->x_outvec;
        x->x_outvec = 0;
    }
    memcpy(outvec, x->x_vec, size);
    outlet_list(x->x_obj.ob_outlet, &s_list, x->x_n, outvec);
    if (reentered)
        t_freebytes(outvec, size);
    else x->x_outvec = outvec;
}

static void pak_pointer(t_pak *x, t_gpointer *gp)
{
    if(x->x_vec->a_type == A_POINTER)
    {
        gpointer_unset(x->x_gpointer);
        *x->x_gpointer = *gp;
        if (gp->gp_stub) gp->gp_stub->gs_refcount++;
        pak_bang(x);
    }
    else pd_error(x, "pak_pointer: wrong type");
}

static void pak_float(t_pak *x, t_float f)
{
    post("float");
    if (x->x_vec->a_type == A_FLOAT)
    {
        x->x_vec->a_w.w_float = f;
        pak_bang(x);
    }
    else pd_error(x, "pak_float: wrong type");
}

static void pak_symbol(t_pak *x, t_symbol *s)
{
    if (x->x_vec->a_type == A_SYMBOL)
    {
        x->x_vec->a_w.w_symbol = s;
        pak_bang(x);
    }
    else pd_error(x, "pak_symbol: wrong type");
}

static void pak_list(t_pak *x, t_symbol *s, int ac, t_atom *av)
{
    obj_list(&x->x_obj, 0, ac, av);
}

static void pak_anything(t_pak *x, t_symbol *s, int ac, t_atom *av)
{
    t_atom *av2 = (t_atom *)getbytes((ac + 1) * sizeof(t_atom));
    int i;
    for (i = 0; i < ac; i++)
        av2[i + 1] = av[i];
    SETSYMBOL(av2, s);
    obj_list(&x->x_obj, 0, ac+1, av2);
    freebytes(av2, (ac + 1) * sizeof(t_atom));
}

static void pak_free(t_pak *x)
{
    t_gpointer *gp;
    int i;
    for(gp = x->x_gpointer, i = x->x_nptr; i--; gp++)
    {
        gpointer_unset(gp);
    }
    freebytes(x->x_vec, x->x_n * sizeof(*x->x_vec));
    freebytes(x->x_outvec, x->x_n * sizeof(*x->x_outvec));
    freebytes(x->x_gpointer, x->x_nptr * sizeof(*x->x_gpointer));
}



extern void pak_setup(void)
{
    t_class* c = NULL;
    
    c = class_new(gensym("spam-inlet"), 0, 0, sizeof(t_pak_inlet), CLASS_PD, 0);
    if(c)
    {
        class_addbang(c,    (t_method)pak_inlet_bang);
        class_addpointer(c, (t_method)pak_inlet_pointer);
        class_addfloat(c,   (t_method)pak_inlet_float);
        class_addsymbol(c,  (t_method)pak_inlet_symbol);
        class_addlist(c,    (t_method)pak_inlet_list);
        class_addanything(c,(t_method)pak_inlet_anything);
    }
    pak_inlet_class = c;
    
    c = class_new(gensym("pak"), (t_newmethod)pak_new, (t_method)pak_free, sizeof(t_pak), CLASS_NOINLET, A_GIMME, 0);
    pak_class = c;
}
