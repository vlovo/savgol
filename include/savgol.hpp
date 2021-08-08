/**
 * C++ Savitzky Golay filter  - Version 0.0.1
 * by M. Leitz MLeitz@boptics.de
 * This Library is licensed under the MIT License
 * https://github.com/vlovo/savgol
 */
#pragma once

#include <array>
#include <algorithm>

namespace filter
{

template <typename InputIt, typename OutputIt, typename Convolution_t>
inline OutputIt savgol(InputIt src_first, InputIt src_last, OutputIt dest_first, Convolution_t smoothOperator)
{

  auto windowBegin = src_first;
  auto windowEnd = src_first + smoothOperator.windowSize();

  auto begin = src_first;

  while (src_first != src_last)
  {

    if (std::distance(begin, src_first) > (smoothOperator.windowSize() / 2 - 1) && std::distance(begin, src_first) < (std::distance(begin, src_last) - (smoothOperator.windowSize() / 2) - 1))
    {
      auto result = accumulate_window(windowBegin, windowEnd, smoothOperator);

      result =  result/smoothOperator.norm();

      *dest_first = result;

      dest_first++;
      src_first++;

      ++windowBegin; // move window forward
      ++windowEnd;   // move window  forward
    }
    else
    {

      *dest_first++ = *src_first++;
    }
  }

  return dest_first;
}

template <typename InputIt, typename ConvOperator>
inline double accumulate_window(InputIt first, InputIt last, ConvOperator op)
{
  double sum=0.0;
  std::size_t i = 0;

  for (; first != last; ++first)
  {
    sum += op(i, *first);
    ++i;
  }

  return sum;
}

// constexpr Map
// inspired by Jason Turner, C++ Weekly Ep 223 on Youtube, my C++ Weekly blog series
// and blog post by  https://xuhuisun.com/post/c++-weekly-2-constexpr-map/

template <typename Key, typename Value, std::size_t Size>
struct c_map
{
  std::array<std::pair<Key, Value>, Size> data;

  [[nodiscard]] constexpr Value at(const Key &key) const
  {
    const auto itr = std::find_if(begin(data), end(data), [&key](const auto &v) { return v.first == key; });
    if (itr != end(data))
    {
      return itr->second;
    }
    else
    {
      throw std::range_error("Not Found");
    }
  }

  [[nodiscard]] constexpr std::size_t count(const Key &key) const
  {
    const auto count = std::count_if(begin(data), end(data), [&key](const auto &v) { return v.first == key; });
    return count;
  }
};

class SmoothQuadCubic
{
public:
  SmoothQuadCubic(int windowSize) : mWindowsize(windowSize)
  {
    if (convolutionWeights.count(windowSize) == 0)
    {
    }
  };

  double operator()(size_t index, double element)
  {
    return static_cast<double>(convolutionWeights.at(mWindowsize)[index + 1] * element);
  };
  auto norm()
  {
    return  (convolutionWeights.at(mWindowsize)[0]);
  };

  int windowSize()
  {
    return mWindowsize;
  }

private:
  int mWindowsize;

  inline static constexpr auto convolutionWeights = c_map<int, std::array<int, 26>, 11>{

      {{{5, {35, -3, 12, 17, 12, -3, 0, 0, 0, 0, 0, 0, 0, 0}},
        {7, {21, -2, 3, 6, 7, 6, 3, -2, 0, 0, 0, 0, 0, 0}},
        {9, {231, -21, 14, 39, 54, 59, 54, 39, 14, -21, 0, 0, 0, 0}},
        {11, {429, -36, 9, 44, 69, 84, 89, 84, 69, 44, 9, -36, 0, 0}},                                                     // Window size 11
        {13, {143, -11, 0, 9, 16, 21, 24, 25, 24, 21, 16, 9, 0, -11}},                                                     // Window size 13
        {15, {1105, -78, -13, 42, 87, 122, 147, 162, 167, 162, 147, 122, 87, 42, -13, -78}},                               // Window size 15
        {17, {323, -21, -6, 7, 18, 27, 34, 39, 42, 43, 42, 39, 34, 27, 18, 7, -6, -21}},                                   // Window size 17
        {19, {2261, -136, -51, 24, 89, 144, 189, 224, 249, 264, 269, 264, 249, 224, 189, 144, 89, 24, -51, -136}},         // Window size 19
        {21, {3059, -171, -76, 9, 84, 149, 204, 249, 284, 309, 324, 329, 324, 309, 284, 249, 204, 149, 84, 9, -76, -171}}, // Window size 21
        {23, {8059, -42, -21, -2, 15, 30, 43, 54, 63, 70, 75, 78, 79, 78, 75, 70, 63, 54, 43, 30, 15, -2, -21, -42}},      // Window size 23
        {25, {5175, -254, -138, -33, 62, 147, 222, 287, 322, 387, 422, 447, 462, 467, 462, 447, 422, 387, 322, 287, 222, 147, 62, -33, -138, -254}}}}

  };

   
};

class SmoothQuarticQuintic
{
public:
  SmoothQuarticQuintic(int windowSize)
      : mWindowsize(windowSize){

        };

  double operator()(size_t index, double element)
  {
    return static_cast<double>(convolutionWeights.at(mWindowsize)[index + 1] * element);
  };
  double norm()
  {
    return static_cast<double>(convolutionWeights.at(mWindowsize)[0]);
  };

  int windowSize()
  {
    return mWindowsize;
  }

private:
  int mWindowsize;

  inline static constexpr auto convolutionWeights = c_map<int, std::array<double, 26>, 10>{

      {{{7, {231, 5, -30, 75, 131, 75, -30, 5, 0, 0, 0, 0, 0, 0}},
        {9, {429, 15, -55, 30, 135, 179, 135, 30, -55, 15, 0, 0, 0, 0}},
        {11, {429, 18, -45, -10, 60, 120, 143, 120, 60, -10, -45, 18, 0, 0}},
        {13, {2431, 110, -198, -160, 110, 390, 600, 677, 600, 390, 110, -160, -198, 110}},
        {15, {46189, 2145, -2860, -2937, -165, 3755, 7500, 10125, 11053, 10125, 7500, 3755, -165, -2937, -2860, 2145}},
        {17, {4199, 195, -195, -260, -117, 135, 415, 660, 825, 883, 825, 660, 415, 135, -117, -260, -195, 195}},
        {19, {7429, 340, -255, -420, -290, 18, 405, 790, 1110, 1320, 1393, 1320, 1110, 790, 405, 18, -290, -420, -255, 340}},
        {21, {260015, 11628, -6460, -13005, -11220, -3940, 6378, 17655, 28190, 36660, 42120, 44003, 42120, 36660, 28190, 17655, 6378, -3940, -11220, -13005, -6460, 11628}},
        {23, {655, 285, -114, -285, -285, -165, 30, 261, 495, 705, 870, 975, 1011, 975, 870, 705, 495, 261, 30, -165, -285, -285, -114, 285}},
        {25, {30015, 1265, -345, -1122, -1255, -915, -255, 590, 1503, 2385, 3155, 3750, 4125, 4253, 4125, 3750, 3155, 2385, 1503, 590, -225, -915, -1255, -1122, -345, 1265}}}}};
};

class DeriveQuadFirst
{
public:
  DeriveQuadFirst(int windowSize)
      : mWindowsize(windowSize){

        };

  double operator()(size_t index, double element)
  {
    return static_cast<double>(convolutionWeights.at(mWindowsize)[index + 1] * element);
  };
  double norm()
  {
    return static_cast<double>(convolutionWeights.at(mWindowsize)[0]);
  };

  int windowSize()
  {
    return mWindowsize;
  }

private:
  int mWindowsize;

  inline static constexpr auto convolutionWeights = c_map<int, std::array<double, 26>, 11>{

      {{{5, {10, 2, 1, 0, -1, -2, 0, 0, 0, 0, 0, 0, 0, 0}},
        {7, {28, 3, 2, 1, 0, -1, -2, -3, 0, 0, 0, 0, 0, 0}},
        {9, {60, 4, 3, 2, 1, 0, -1, -2, -3, -4, 0, 0, 0, 0}},
        {11, {110, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, 0, 0, 0, 0}},
        {13, {182, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, 0, 0}},
        {15, {280, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7}},
        {17, {408, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8}},
        {19, {570, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9}},
        {21, {770, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10}},
        {23, {1012, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11}}, 
        {25, {1300, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12}}}}

  };
};

} // namespace filter
