/*
 // Copyright (c) 2016 Pierre Guillot.
 // For information on usage and redistribution, and for a DISCLAIMER OF ALL
 // WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include <m_pd.h>

static t_class *rotilt_class;

typedef struct _rotilt
{
    t_object f_obj;
    t_float  f_f;
} t_rotilt;

static void *rotilt_new()
{
    return pd_new(rotilt_class);
}
/*
t_int *rotilt_perform(t_int *w)
{
    int i;
    t_sample  *inw = (t_sample *)(w[1]);
    t_sample  *inx = (t_sample *)(w[2]);
    t_sample  *iny = (t_sample *)(w[3]);
    t_sample  *inz = (t_sample *)(w[4]);
    
    t_sample  *inax = (t_sample *)(w[5]);
    t_sample  *inay = (t_sample *)(w[6]);
    t_sample  *inaz = (t_sample *)(w[7]);
    
    t_sample  *outw = (t_sample *)(w[8]);
    t_sample  *outx = (t_sample *)(w[9]);
    t_sample  *outy = (t_sample *)(w[10]);
    t_sample  *outz = (t_sample *)(w[11]);
    int          n = (int)(w[12]);
    
    while(n--)
    {
        (*outx++) = sample2 * ((cosf(sample10)-(sample3 * (sinf(sample10))))) * ((sample2*cosf(sample11))-(sample4*(sinf(sample11)));
    }
    
                                                                               
    (*APout1++) = sample1 * 1; //W
    (*APout2++) = sample2 * ((cosf(sample10)-(sample3 * (sinf(sample10))))) * ((sample2*cosf(sample11))-(sample4*(sinf(sample11)));
                                                                               (*APout3++) = sample2 * ((sinf(sample10)+(sample3 * (cosf(sample10))))) * ((sample3*cosf(sample12))-(sample4*(sinf(sample12)));
                                                                                                                                                          (*APout4++) = sample2 * ((sinf(sample11)+(sample4 * (cosf(sample11))))) * ((sample3*sinf(sample12))+(sample4*(cosf(sample12)));
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
    
    return (w+13);
}

 */

void rotilt_dsp(t_rotilt *x, t_signal **sp)
{
    /*
    dsp_add(rotilt_perform, 12, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, // Input Harmonics
            sp[4]->s_vec, sp[5]->s_vec, sp[6]->s_vec, // Inputs Angles
            sp[7]->s_vec, sp[8]->s_vec, sp[9]->s_vec, sp[10]->s_vec, // Output Harmonics
            sp[0]->s_n);
     */
}


extern void rotilt_tilde_setup(void)
{
    t_class* c = class_new(gensym("rotilt~"), (t_newmethod)rotilt_new, NULL, sizeof(t_rotilt), CLASS_DEFAULT, 0);
    if(c)
    {
        CLASS_MAINSIGNALIN(c, t_rotilt, f_f);
        class_addmethod(c, (t_method)rotilt_dsp, gensym("dsp"), A_CANT);
        
    }
    rotilt_class = c;
}
