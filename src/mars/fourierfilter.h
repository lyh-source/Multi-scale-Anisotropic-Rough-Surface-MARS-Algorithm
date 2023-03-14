#pragma once
#include<iostream>
#include<opencv2/opencv.hpp>
using namespace std;
class FourierFilter {
public:
	FourierFilter() {};
	~FourierFilter() {};
	void setValue(const int& option_, const int& cutoff_);
	int option, cutoff;


    void fft2(const cv::Mat& src, cv::Mat& Fourier)
    {
        int mat_type = src.type();
        assert(mat_type < 15); //Unsupported Mat datatype

        if (mat_type < 7)
        {
            cv::Mat planes[] = { cv::Mat_<double>(src), cv::Mat::zeros(src.size(),CV_64F) };
            merge(planes, 2, Fourier);
            cv::dft(Fourier, Fourier);
        }
        else // 7<mat_type<15
        {
            cv::Mat tmp;
            dft(src, tmp);
            vector<cv::Mat> planes;
            split(tmp, planes);
            magnitude(planes[0], planes[1], planes[0]); //Change complex to magnitude
            Fourier = planes[0];
        }
    }

    void ifft2(const cv::Mat& src, cv::Mat& Fourier)
    {
        int mat_type = src.type();
        assert(mat_type < 15); //Unsupported Mat datatype

        if (mat_type < 7)
        {
            cv::Mat planes[] = { cv::Mat_<double>(src), cv::Mat::zeros(src.size(),CV_64F) };
            merge(planes, 2, Fourier);
            dft(Fourier, Fourier, cv::DFT_INVERSE + cv::DFT_SCALE, 0);
        }
        else // 7<mat_type<15
        {
            cv::Mat tmp;
            dft(src, tmp, cv::DFT_INVERSE + cv::DFT_SCALE, 0);
            //vector<Mat> planes;
            //split(tmp, planes);

            //magnitude(planes[0], planes[1], planes[0]); //Change complex to magnitude
            //Fourier = planes[0];
            Fourier = tmp;
        }
    }

};