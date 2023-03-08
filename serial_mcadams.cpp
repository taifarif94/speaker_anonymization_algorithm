#include <iostream>
#include "Audiofile.h"
#include <cmath>
#include <string>
#include <complex.h>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <algorithm>

using std::cout;
using std::endl;
using namespace std;

#define FILENAME "newlpc.txt"
#define FILE_FRAME_REC "frame_rec.txt"

float *hanning(int N, short itype);

//%Parameters%
const uint8_t winLengthinms = 20;
const uint8_t shiftLengthinms = 10;
const uint8_t lp_order = 20;
const float mcadams = 0.8;



int main()

{
    // File Parameters
    string lpc_coefficients;
    /* Open the file for reading */
    char *line_buf = NULL;
    size_t line_buf_size = 0;
    int line_count = 0;
    ssize_t line_size;
    FILE *fp = fopen(FILENAME, "r");
    if (!fp)
    {
        fprintf(stderr, "Error opening file '%s'\n", FILENAME);
        return EXIT_FAILURE;
    }

    /* Get the first line of the file. */
    line_size = getline(&line_buf, &line_buf_size, fp);

    std::complex<double> iota(0.0, 1.0);
    std::complex<double> iotan(0.0, -1.0);

    std::cout << std::fixed;

    //Audio file import settings

    AudioFile<double> audioFile;
    AudioFile<long double> savedFile;
    savedFile.setNumChannels(1);

    audioFile.load("in.wav");
    uint8_t channel = 0;
    uint32_t sig = audioFile.getNumSamplesPerChannel();

   savedFile.setNumSamplesPerChannel(sig);

    uint32_t fs = audioFile.getSampleRate();

    for (int i = 0; i < sig; i++)
    {
        audioFile.samples[channel][i] = audioFile.samples[channel][i] + std::numeric_limits<double>::epsilon();
    }

    //simulation parameters

    uint16_t winlen = floor(winLengthinms * 0.001 * fs);
    cout << "winlen: " << winlen << endl;
    uint16_t shift = floor(shiftLengthinms * 0.001 * fs);
    cout << "shift: " << shift << endl;
    float length_sig = sig;
    uint16_t length_sigs = sig;
    cout << "length_sig: " << length_sig << endl;
    float duration_sig = length_sig / fs;
    cout << "duration_sig: " << duration_sig << endl;

    //FFT processing parameters

    uint16_t NFFT = pow(2, ((uint16_t)((ceil(log2(winlen))))));
    cout << "NFFT: " << NFFT << endl;
    float *wPR = hanning(winlen, 0);
    cout << "wPR[319]: " << wPR[319] << endl;
    // This array starts with 0 and is equivalent to position 'n+2' of Matlab in the start but equivalent for later??
    // it goes from '0' to matlab-1. Matlab goes from '1' to matlab.
    double sum_wPR = 0;

    for (int i = 0; i < winlen; i++)
    {
        sum_wPR = sum_wPR + *(wPR + i);
    }

    double K = sum_wPR / shift;
    // A slight difference of K between Matlab and this 'K'
    cout << "K: " << K << endl;

    double *win;
    win = new double[winlen];

    for (int i = 0; i < winlen; i++)
    {
        win[i] = sqrt((*(wPR + i) / K));
    }
    cout << "win[319]: " << win[319] << endl;

    uint16_t Nframes = 1 + floor((length_sig - winlen) / shift);
    double *sig_rec;
    sig_rec = new double[length_sigs];

    for (long i = 0; i < length_sigs; i++)
    {
        sig_rec[i] = 0;
    }

    cout << "Nframes: " << Nframes << endl;

    std::vector<long double> _frame_rec_super(length_sigs);
    for (int x = 0; x < length_sigs; x++)
        _frame_rec_super[x] = 0;

    for (int m = 1; m < Nframes; m++)
    {
        uint16_t *index;
        uint16_t minimum = m * shift + winlen <= length_sigs ? m * shift + winlen : length_sigs;
        uint16_t minimum_index = minimum - m * shift;
        index = new uint16_t[minimum];
        // error for the index and frame:
        for (int j = 0; j < minimum_index; j++)

        {
            index[j] = ((m * shift) + j + 1);
        }

        double *frame;
        frame = new double[minimum_index];

        for (int k = 0; k < winlen; k++)
        {
            frame[k] = audioFile.samples[channel][k] * win[k];
        }

        while (line_size >= 0)
        {
            
            /* Increment our line count */
            line_count++;
            /* Get the next line */
            line_size = getline(&line_buf, &line_buf_size, fp);
            if (line_count == 1)
            {
                lpc_coefficients = line_buf;

            }

        }
        std::stringstream ss(lpc_coefficients);
        std::vector<string> vect;
        std::vector<long double> lpc_coeff;
        while (ss.good())
        {
            string substr;
            getline(ss, substr, ',');
            vect.push_back(substr);

            if (!(substr.empty()))
                lpc_coeff.push_back(std::stold(substr));
        }


        /* Free the allocated line buffer */
        free(line_buf);
        line_buf = NULL;

        // code for reading and saving file contents of file "FRAME_REC.txt"

        string _frame_rec;
        /* Open the file for reading */
        char *line_buf_frame_rec = NULL;
        size_t line_buf_size_frame_rec = 0;
        int line_count_frame_rec = 0;
        ssize_t line_size_frame_rec;
        FILE *fp_frame_rec = fopen(FILE_FRAME_REC, "r");
        if (!fp_frame_rec)
        {
            fprintf(stderr, "Error opening file '%s'\n", FILENAME);
            return EXIT_FAILURE;
        }

        /* Get the first line of the file. */
        line_size_frame_rec = getline(&line_buf_frame_rec, &line_buf_size_frame_rec, fp_frame_rec);


        /* Loop through until we are done with the file. */
        while (line_size_frame_rec >= 0)
        {
            /* Increment our line count */
            ++line_count_frame_rec;
            if (line_count_frame_rec == m)
            {
                _frame_rec = line_buf_frame_rec;

            }
            line_size_frame_rec = getline(&line_buf_frame_rec, &line_buf_size_frame_rec, fp_frame_rec);
        }

        std::stringstream ss_frame_rec(_frame_rec);
        std::vector<string> vect_frame_rec;
        std::vector<long double> _frame_rec_d;

        while (ss_frame_rec.good())
        {
            string substr_frame_rec;
            getline(ss_frame_rec, substr_frame_rec, ',');
            vect_frame_rec.push_back(substr_frame_rec);
            if (!(substr_frame_rec.empty()))
                _frame_rec_d.push_back(std::stold(substr_frame_rec));
        }

        for (std::size_t i = 0; i < _frame_rec_d.size(); i++)
        {
            _frame_rec_d[i] = _frame_rec_d[i] * win[i];
        }

        /* Free the allocated line buffer */
        free(line_buf_frame_rec);
        line_buf_frame_rec = NULL;

        /* Close the file now that we are done with it */
        fclose(fp_frame_rec);

        std::vector<long> outindex;
        outindex.clear();

        for (int p = (m * shift); p < ((m * shift) + _frame_rec_d.size()); p++)
        {
            outindex.push_back(p);
        }

        for (int e = 0; e < outindex.size(); e++)
        {
            _frame_rec_super[outindex[e]] = _frame_rec_super[outindex[e]] + _frame_rec_d[e];

            if(m == 125 && e <=320 && e>=311)
            std::cout<<_frame_rec_super[outindex[e]]<<"  ";

        }

    }

    std::cout << _frame_rec_super[30235] << std::endl;

    long double smallest_element = _frame_rec_super[0];       // let, first element is the smallest one
    long double largest_element = _frame_rec_super[0];        // also let, first element is the biggest one


     for (int y = 0; y < _frame_rec_super.size(); y++) 
    {
        if (_frame_rec_super[y] < 0)
        _frame_rec_super[y] = (_frame_rec_super[y]*-1);
   
    }

    for (int i = 1; i < _frame_rec_super.size(); i++) 
    {
        if (_frame_rec_super[i] < smallest_element)
        {
            smallest_element = _frame_rec_super[i];
        }
        if (_frame_rec_super[i] > largest_element)
        {
            largest_element = _frame_rec_super[i];
        }
    }

    for (int u = 0; u < _frame_rec_super.size(); u++) 
    {
        
        _frame_rec_super[u] = (_frame_rec_super[u]/largest_element);
   
    }

    std::cout<<"Largest Value: "<<largest_element<<std::endl;
    std::cout<<"Smallest Value: "<<smallest_element<<std::endl;


    for (int i = 0; i < audioFile.getNumSamplesPerChannel(); i++)
        {
            for (int channel = 0; channel < audioFile.getNumChannels(); channel++)
            {
                
                audioFile.samples[channel][i] = _frame_rec_super[i];
            }
        }

    audioFile.save ("sample_output.wav", AudioFileFormat::Wave);

    std::cout<<"Sample Rate: "<<audioFile.getSampleRate()<<std::endl;

    // complex array comntains dummy poles
    std::complex<double> complex_array[100];

    for (int i = 0; i < 100; i++)
    {
        if (i < 50)
            complex_array[i] = std::complex<double>(0.0, 0.0);
        else
            complex_array[i] = std::complex<double>(50.0, 50.0);
    }

    std::vector<int> ind_imag;
    std::vector<int> ind_imag_con;

    // size of i for the ind_imag for loop must be determined at runtime
    for (int i = 0; i < 100; i++)
    {
        if (imag(complex_array[i]) != 0)
        {
            ind_imag.push_back(i);

        }
    }

    for (int i = 0; i < ind_imag.size(); i++)
    {
        if (i == 0)
        {
            ind_imag_con.push_back(ind_imag[i]);
        }
        if (i % 2 != 0)
        {
            ind_imag_con.push_back(ind_imag[i]);
        }
    }

    double new_angles_check;
    std::vector<double> new_angles;
    std::complex<double> complex_for_phase;

    complex_for_phase.real(1.7552);
    complex_for_phase.imag(0.9589);

    for (int i = 0; i < ind_imag_con.size(); i++)

    {
        new_angles_check = pow(arg(complex_array[ind_imag_con[i]]), mcadams);
        if (new_angles_check > 3.141592654)
            new_angles.push_back(3.141592654);
        else if (new_angles_check < 0)
            new_angles.push_back(0);
        else
            new_angles.push_back(new_angles_check);
    }

    std::complex<double> new_poles[100];

    for (int i = 0; i < 100; i++)
    {
        new_poles[i] = complex_array[i];
    }

    for (int k = 0; k < ind_imag_con.size() - 1; k++)

    {
        new_poles[ind_imag_con[k]] = abs(complex_array[ind_imag_con[k]]) * exp(iota * new_angles[k]);
        new_poles[ind_imag_con[k + 1]] = abs(complex_array[ind_imag_con[k + 1]]) * exp(iotan * new_angles[k]);
    }


    fclose(fp);


    return 0;
}

/*  function w = hanning(varargin)
%   HANNING   Hanning window.
%   HANNING(N) returns the N-point symmetric Hanning window in a column
%   vector.  Note that the first and last zero-weighted window samples
%   are not included.
%
%   HANNING(N,'symmetric') returns the same result as HANNING(N).
%
%   HANNING(N,'periodic') returns the N-point periodic Hanning window,
%   and includes the first zero-weighted window sample.
%
%   NOTE: Use the HANN function to get a Hanning window which has the
%          first and last zero-weighted samples.ep
    itype = 1 --> periodic
    itype = 0 --> symmetric
    default itype=0 (symmetric)

    Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.3 $  $Date: 2007/12/14 15:05:04 $
*/

float *hanning(int N, short itype)
{
    int half, i, idx, n;
    float *w;

    w = (float *)calloc(N, sizeof(float));
    memset(w, 0, N * sizeof(float));

    if (itype == 1) // periodic function
        n = N - 1;
    else
        n = N;

    if (n % 2 == 0)
    {
        half = n / 2;
        for (i = 0; i < half; i++) // CALC_HANNING   Calculates Hanning window samples.
            w[i] = 0.5 * (1 - cos(2 * M_PI * (i + 1) / (n + 1)));

        idx = half - 1;
        for (i = half; i < n; i++)
        {
            w[i] = w[idx];
            idx--;
        }
    }
    else
    {
        half = (n + 1) / 2;
        for (i = 0; i < half; i++) // CALC_HANNING   Calculates Hanning window samples.
            w[i] = 0.5 * (1 - cos(2 * M_PI * (i + 1) / (n + 1)));

        idx = half - 2;
        for (i = half; i < n; i++)
        {
            w[i] = w[idx];
            idx--;
        }
    }

    if (itype == 1) // periodic function
    {
        for (i = N - 1; i >= 1; i--)
            w[i] = w[i - 1];
        w[0] = 0.0;
    }
    return (w);
}

