#include <opencv2/core.hpp>

#include "highmap.hpp"

int main(void)
{
  cv::Mat     mat = cv::Mat::zeros(256, 256, CV_32FC1);
  hmap::Array z = hmap::cv_mat_to_array(mat);
  return 0;
}
