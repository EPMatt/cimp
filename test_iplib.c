#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
#include "time.h"

#define ROUNDS_NUM 10 /* number of tests iterations */
#define MAT_W 100
#define MAT_H 200
#define MAT_CHANNELS 3
#define TOTAL_TESTS_NUM MAT_CHANNELS + 23

int main()
{
    int round;
    ip_mat **round_arr;
    /* init random */
    srand(time(NULL));
    for (round = 1; round <= ROUNDS_NUM; round++)
    {
        channel_t a, b;
        ip_mat **pt, *original_pt, *a_pt, *b_pt, *a_filter, *padded;
        int row, col, ch, size_w, size_h;
        float foo_mean;
        printf("TEST No.%d: Start...\n", round);
        /* generate random sizes */
        size_h = rand() % MAT_H + 1;
        size_w = rand() % MAT_W + 1;
        printf("\tGenerated random sizes ->  height: %d, width: %d\n", size_h, size_w);
        round_arr = (ip_mat **)malloc(sizeof(ip_mat *) * (TOTAL_TESTS_NUM));
        pt = round_arr;
        /* execute all the exposed library functions (those declared in the header file) */
        /* PART 1: mem */
        printf("\tip_mat_create...\n");
        *pt = ip_mat_create(size_h, 100, MAT_CHANNELS, 0.0);
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
        original_pt = *pt;
        pt++;
        printf("\tip_mat_init_random...\n");
        ip_mat_init_random(original_pt, 10.0, 0.5);
        printf("\tip_mat_copy...\n");
        *pt = ip_mat_copy(original_pt);
        pt++;
        printf("\tip_mat_subset...\n");
        *pt = ip_mat_subset(original_pt, 0, original_pt->h, 0, original_pt->w);
        a_pt = original_pt;
        b_pt = original_pt;
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
        *pt = ip_mat_mul_scalar(original_pt, 10.0);
        pt++;
        printf("\tip_mat_add_scalar...\n");
        *pt = ip_mat_add_scalar(original_pt, 10.0);
        pt++;
        printf("\tip_mat_mean...\n");
        *pt = ip_mat_mean(a_pt, b_pt);
        pt++;
        /* PART 2 */
        printf("\tip_mat_to_gray_scale...\n");
        *pt = ip_mat_to_gray_scale(original_pt);
        pt++;
        printf("\tip_mat_blend...\n");
        *pt = ip_mat_blend(a_pt, b_pt, 0.5);
        pt++;
        printf("\tip_mat_brighten...\n");
        *pt = ip_mat_brighten(original_pt, 0.5);
        pt++;
        printf("\tip_mat_corrupt...\n");
        *pt = ip_mat_corrupt(original_pt, 2.0);
        pt++;
        /* PART 3 */
        printf("\tip_mat_padding...\n");
        *pt = ip_mat_padding(original_pt, 1, 1);
        padded = *pt;
        pt++;
        /* first make some filters */
        printf("\tcreate_sharpen_filter...\n");
        *pt = create_sharpen_filter();
        a_filter = *pt;
        pt++;
        printf("\tcreate_edge_filter...\n");
        *pt = create_edge_filter();
        pt++;
        printf("\tcreate_emboss_filter...\n");
        *pt = create_emboss_filter();
        pt++;
        printf("\tcreate_average_filter...\n");
        *pt = create_average_filter(3, 3, 3);
        pt++;
        printf("\tcreate_gaussian_filter...\n");
        *pt = create_gaussian_filter(3, 3, 3, 0.5);
        pt++;
        /* test convolve: all filters */
        printf("\tip_mat_convolve (sharpen_filter)...\n");
        *pt = ip_mat_convolve(padded, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (edge_filter)...\n");
        *pt = ip_mat_convolve(padded, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (emboss_filter)...\n");
        *pt = ip_mat_convolve(padded, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (average_filter)...\n");
        *pt = ip_mat_convolve(padded, *(pt - 5));
        pt++;
        printf("\tip_mat_convolve (gaussian_filter)...\n");
        *pt = ip_mat_convolve(padded, *(pt - 5));
        /* clamp & rescale on last ip_mat */
        printf("\tclamp...\n");
        clamp(*pt, 50.0, 75.0);
        printf("\trescale...\n");
        rescale(*pt, 100.0);
        /* HELPERS */
        printf("\tget_channel...\n");
        a = get_channel(padded, 0);
        b = get_channel(padded, 1);
        printf("\tchannel_puts...\n");
        channel_puts(b, a, 0, 0);
        ip_mat_puts(b_pt, a_pt, 0, 0, NO_COMPUTE_STATS);
        printf("\tnot_null_ip_mat...\n");
        not_null_ip_mat(a_pt);
        printf("\tmin...\n");
        min(a_pt, 0);
        printf("\tmax...\n");
        max(a_pt, 0);
        printf("\tmean...\n");
        foo_mean = mean(a_pt, 0);
        printf("\trestrict_val...\n");
        restrict_val(foo_mean, MIN_PIXEL_FLOAT, MAX_PIXEL_FLOAT);
        printf("\tmean_pixel_channel...\n");
        mean_pixel_channel(original_pt, 0, 0);
        printf("\tequal_dimension...\n");
        equal_dimension(a_pt, b_pt);
        printf("\tconvolve_channel...\n");
        convolve_channel(a, get_channel(a_filter, 0), 0, 0);
        /* free everything */
        printf("\tip_mat_free...\n");
        for (pt = round_arr; pt < round_arr + TOTAL_TESTS_NUM; pt++)
            ip_mat_free(*pt);
        printf("\tall tests done! free the round array...\n");
        free(round_arr);
        printf("TEST No.%d: Done!\n", round);
    }
}