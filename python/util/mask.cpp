    // #include <cstdint>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

# define M_PI           3.14159265358979323846

typedef enum {
    CHIRP = 0,
    SPOTS = 1,
    BARS = 2,
    NONE = 255
} trial_t;


void chirp(uint8_t* pattern, double intensity, uint64_t frame_rate) {
    // goes up to 2098...
    printf("frame_rate is: %lld\n", frame_rate);

    double iintensity = 255 * intensity; // hmm....

    memset(&pattern[0 * frame_rate], 0,                             2 * frame_rate); //pre
    memset(&pattern[2 * frame_rate], (uint8_t) round(iintensity),   3 * frame_rate); // pos
    memset(&pattern[5 * frame_rate], 0,                             3 * frame_rate); // neg
    memset(&pattern[8 * frame_rate], (uint8_t) round(iintensity/2), 2 * frame_rate); // inter    
    
    double freq = 0.0;
    for (auto i = 0; i < 8*frame_rate; i++) {
        freq += 8.0/frame_rate * i / (8.0*frame_rate - 1);
        pattern[10 * frame_rate + i] = (uint8_t) round((-sin(2 * M_PI * freq + M_PI)/2 + .5) * iintensity); //freq chirp
    }

    
    memset(&pattern[18 * frame_rate], (uint8_t) round(iintensity/2),2 * frame_rate); // inter

    // double amp = 0;
    double phase = 0; //1.0 / frame_rate;
    for (auto i = 0; i < 8*frame_rate; i++) {
        phase += 1.0 / frame_rate;
        pattern[20 * frame_rate + i] = (uint8_t) round((-i / (8.0*frame_rate - 1) * sin(4 * M_PI * phase + M_PI)/2  + .5) * iintensity); //amp chirp
    }

    //
    memset(&pattern[28 * frame_rate], (uint8_t) round(iintensity/2),2 * frame_rate); // inter
    memset(&pattern[30 * frame_rate], 0,                            2098 - 30 * frame_rate); // tail
}

void mask(
    double* y,
    double* x,
    double* w,

    int width,
    int n_trials,

    uint64_t* frames,
    uint64_t* n_frames,
    uint64_t* flips,
    uint64_t* trial_type,
    uint64_t* pre_frames,
    uint64_t* stim_frames,
    uint64_t* tail_frames,

    double* frame_rate,
    double* intensity
    // uint8_t* intensity
) {

    uint8_t pattern[2100]; //constant is a bit weird... TODO: fix this
    //TODO: it would make more sense / be a bit faster to have a different pattern for each trial type?
    // just switch between them using a pointer...
    // update as necessary


    uint64_t k = 0;
    uint64_t l = 0;

    double last_intensity = 0.0f;
    uint64_t last_trial = NONE;
    uint64_t last_pre = 0;
    uint64_t last_stim = 0;  
    uint64_t last_tail = 0;  

    printf("width is: %d\n", width);

    printf("%d trials\n", n_trials);
    for (auto i=0; i<n_trials; i++) { // for trial in trials
        if ((i % 25) == 0) printf(".");
        if ((trial_type[i] != last_trial) || (intensity[i] != last_intensity) || ((trial_type[i] == SPOTS) && ((pre_frames[i] != last_pre) || (stim_frames[i] != last_stim) || (tail_frames[i] != last_tail)) )) {
            switch (trial_type[i]){
                    case CHIRP:
                        printf("generating chirp pattern with frame rate = %lf, intensity = %lf\n", frame_rate[i], intensity[i]);
                        chirp(&pattern[0], intensity[i], (uint64_t) frame_rate[i]);
                        break;
                    case SPOTS:
                        {
                        //only does the first one??
                        auto ii = 0;
                        do {
                            memset(&pattern[ii],0,pre_frames[i]);
                            ii += pre_frames[i];
                            memset(&pattern[ii], (uint8_t) round(intensity[i] * 255), stim_frames[i]);
                            ii += stim_frames[i];
                            memset(&pattern[ii], 0,tail_frames[i]);
                            ii += tail_frames[i];                            

                        } while (ii < 2100);

                        // memset(&pattern[0], 0, pre_frames[i]);
                        // memset(&pattern[pre_frames[i]], (uint8_t) round(intensity[i] * 255), stim_frames[i]);
                        // memset(&pattern[pre_frames[i] + stim_frames[i]], 0, 2098 - pre_frames[i] - stim_frames[i]);   // ??
                        
                        last_pre = pre_frames[i];
                        last_stim = stim_frames[i];    
                        last_tail = tail_frames[i];
                        }                    
                        break;        
                                 
                    case BARS:
                        {
                        // memset(&pattern[0], 0, 15);
                        // memset(&pattern[15], (uint8_t) round(intensity[i] * 255), 180);


                        // pre_i = 0
                        // on_i = 15
                        // off_i = 195
                        // tail_i = 210
                        auto ii = 0;
                        do {
                            memset(&pattern[ii],0, 15);
                            ii += 15;
                            memset(&pattern[ii], (uint8_t) round(intensity[i] * 255), 180);
                            ii += 180;
                            memset(&pattern[ii], 0, 15);
                            ii += 15;                            

                        } while (ii < 2100);
                        }
                        break;

            }
            last_trial = trial_type[i];
            last_intensity = intensity[i];

        }
        uint64_t a = flips[l + i]; // add i, because there's one more flip per trial than frame
        uint64_t lasta = a;
        for (auto j=1; j <= n_frames[i]; j++) { // for frame in trial

            if (frames[l + j - 1]>=2098) break;

            uint64_t b = flips[l+i+j]; // 

            bool debugPrint = false;
            
            // if (((((a-1) % width)==0) && (j>1) && (pattern[frames[l+j-2]]!=0) && (pattern[frames[l+j-1]]!=0)) 
            //     || ((((b-1) % width)==0) && (j<n_frames[i]) && (pattern[frames[l + j]]!=0) && (pattern[frames[l+j-1]]!=0))) {
            if ((((a-1) % width)==0) || (((b-1) % width)==0)) {
                // if either flip is on a boundary, set to nan (conservative?)
                // however, we only do this if either frame is nonzero

                // NOTE: in ieee, all ones is nan, so we can use memset
                memset(&y[a], 255, (b - a)*8); //number of bytes!!!

                if (debugPrint) printf("!");
                goto next;
            }

            double* wi = &w[pattern[frames[l + j - 1]] * 256];

            double x0 = x[a];
            double xn = 256.0 / (x[b] - x[a]);

            bool ind1 = false;
            bool ind8 = false;

            if ((a >= 3428*32*256) && (a<=3488*32*256)) { //first imaging frame of light stimulus for chirp
                debugPrint = true;
                printf("a_f: %llu, a_x: %llu, a_y: %i\n", a/(32*256), a%256, (a/256)%32);
                printf("\tb_f: %llu, b_x: %llu, b_y: %i\n", b/(32*256), b%256, (b/256)%32);
                printf("\tp: %i\n", pattern[frames[l + j - 1]]);

            }

            // NOTE: should be a lot faster to use blas??
            for (uint64_t  m=a; m<b; m++) {
                uint8_t ind = (uint8_t) ((x[m] - x0) * xn);
                y[m] *= wi[ind];

                ind1 |= (ind==1);
                ind8 |= (ind==8);
            }

            // if ((!ind1) && (!ind8) && (j>1) && (pattern[frames[l+j-2]]!=0) && (pattern[frames[l+j-1]]!=0)) {
            if ((!ind1) && (!ind8)) {
                //for some reason, we have an alignment issue with this case
                //
                memset(&y[lasta], 255, (b-lasta)*8);                
                if (debugPrint) printf("?");
            } else if (debugPrint) printf(".");
next:
            a = b;
            lasta = a;
        }
        l+=n_frames[i];
    }
    printf("\n");

}