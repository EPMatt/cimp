#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

#define ROUNDS_NUM 1 /* number of tests iterations */
#define MAT_W 10
#define MAT_H 20
#define MAT_CHANNELS 3

int main()
{
    int round;
    ip_mat **round_arr;
    for (round = 0; round < ROUNDS_NUM; round++)
    {
        ip_mat **pt, *a_pt, *b_pt;
        int row, col, ch;
        round_arr = (ip_mat **)malloc(sizeof(ip_mat *) * (18 + MAT_CHANNELS));
        pt = round_arr;
        printf("TEST No.%d: Start...\n", round);
        /* execute all the exposed library functions (those declared in the header file) */
        /* PART 1: mem */
        printf("\tip_mat_create...\n");
        *pt = ip_mat_create(MAT_H, MAT_W, MAT_CHANNELS, 0.0);
        /* test get_val & set_val */
        printf("\tget_val & set_val...\n");
        for (ch = 0; ch < (*pt)->k; ch++)
        {
            for (row = 0; row < (*pt)->h; row++)
            {
                for (col = 0; col < (*pt)->w; col++)
                {
                    /* silly get & set test */
                    set_val(*pt, row, col, ch, 1.0);
                }
            }
        }
        printf("\tcompute_stats...\n");
        compute_stats(*pt);
        printf("\tip_mat_init_random...\n");
        ip_mat_init_random(*pt, 10.0, 0.5);
        pt++;
        printf("\tip_mat_copy...\n");
        *pt = ip_mat_copy(*(pt - 1));
        pt++;
        printf("\tip_mat_subset...\n");
        *pt = ip_mat_subset(*(pt - 1), 0, MAT_H, 0, MAT_W);
        pt++;
        a_pt = *(pt - 3);
        b_pt = *(pt - 2);
        pt++;
        /* test concat on all channels */
        printf("\tip_mat_concat...\n");
        for (ch = 0; ch < a_pt->k; ch++)
        {
            *pt = ip_mat_concat(a_pt, b_pt, ch);
            pt++;
        }
        /* PART 1: Math */
        printf("\tip_mat_sum...\n");
        *pt = ip_mat_sum(a_pt, b_pt);
        pt++;
        printf("\tip_mat_sub...\n");
        *pt = ip_mat_sub(a_pt, b_pt);
        pt++;
        printf("\tip_mat_mul_scalar...\n");
        *pt = ip_mat_mul_scalar(*(pt - 1), 10.0);
        pt++;
        printf("\tip_mat_add_scalar...\n");
        *pt = ip_mat_add_scalar(*(pt - 1), 10.0);
        pt++;
        printf("\tip_mat_mean...\n");
        *pt = ip_mat_mean(a_pt, b_pt);
        pt++;
        /* PART 2 */
        printf("\tip_mat_to_gray_scale...\n");
        *pt = ip_mat_to_gray_scale(*(pt - 1));
        pt++;
        printf("\tip_mat_blend...\n");
        *pt = ip_mat_blend(a_pt, b_pt, 0.5);
        pt++;
        printf("\tip_mat_brighten...\n");
        *pt = ip_mat_brighten(*(pt - 1), 0.5);
        pt++;
        printf("\tip_mat_corrupt...\n");
        *pt = ip_mat_corrupt(*(pt - 1), 2.0);
        /* PART 3 */
        /* first make some filters */
        printf("\tcreate_sharpen_filter...\n");
        *pt = create_sharpen_filter();
        pt++;
        printf("\tcreate_edge_filter...\n");
        *pt = create_edge_filter();
        pt++;
        printf("\tcreate_emboss_filter...\n");
        *pt = create_emboss_filter();
        pt++;
        printf("\tcreate_average_filter...\n");
        *pt = create_average_filter(5, 5, 3);
        pt++;
        printf("\tcreate_gaussian_filter...\n");
        *pt = create_gaussian_filter(5, 5, 3, 0.5);
        pt++;
        /* test convolve: all filters */
        printf("\tip_mat_convolve (sharpen_filter)...\n");
        *pt = ip_mat_convolve(a_pt, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (edge_filter)...\n");
        *pt = ip_mat_convolve(a_pt, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (emboss_filter)...\n");
        *pt = ip_mat_convolve(a_pt, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (average_filter)...\n");
        *pt = ip_mat_convolve(a_pt, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (gaussian_filter)...\n");
        *pt = ip_mat_convolve(a_pt, *(pt - 5));
        pt++;
        /* clamp & rescale on last ip_mat */
        printf("\trescale...\n");
        rescale(*pt, 100.0);
        printf("\tclamp...\n");
        clamp(*pt, 50.0, 75.0);
        /* free everything */
        printf("\tip_mat_free...\n");
        for (pt = round_arr; pt < round_arr + (18 + MAT_CHANNELS); pt++)
            ip_mat_free(*pt);
        printf("\tall tests done! free the round array...\n");
        free(round_arr);
        printf("TEST No.%d: Done!\n", round);
    }
}