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

  if (src_first == src_last) return dest_first;

  auto windowBegin = src_first;
  auto windowEnd = src_first + smoothOperator.windowSize();

  const auto data_begin = src_first;

  const auto border_size = smoothOperator.windowSize() / 2;

  while (src_first != src_last)
  {
    if (std::distance(src_first, src_last) <= border_size || std::distance(data_begin, src_first) < border_size)
    {
      *dest_first = *src_first;
    }
    else
    {

      auto result = accumulate_window(windowBegin, windowEnd, smoothOperator);

      result = result / smoothOperator.norm();

      *dest_first = result;

      if (windowEnd != src_last)
      {
        ++windowBegin; // move window forward
        ++windowEnd;   // move window  forward
      }
    }

    ++dest_first;
    ++src_first;
  }

  return dest_first;
}

template <typename InputIt, typename ConvOperator>
inline double accumulate_window(InputIt first, InputIt last, ConvOperator op)
{
  double sum = 0.0;
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
      throw std::invalid_argument("window size not supported");
    }
  };

  double operator()(size_t index, double element)
  {
    return static_cast<double>(convolutionWeights.at(mWindowsize)[index + 1] * element);
  };
  auto norm()
  {
    return (convolutionWeights.at(mWindowsize)[0]);
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
        {23, {805, -42, -21, -2, 15, 30, 43, 54, 63, 70, 75, 78, 79, 78, 75, 70, 63, 54, 43, 30, 15, -2, -21, -42}},       // Window size 23
        {25, {5175, -254, -138, -33, 62, 147, 222, 287, 343, 387, 422, 447, 462, 467, 462, 447, 422, 387, 343, 287, 222, 147, 62, -33, -138, -254}}}}

  };
};

class SmoothQuarticQuintic
{
public:
  SmoothQuarticQuintic(int windowSize) : mWindowsize(windowSize)
  {
    if (convolutionWeights.count(windowSize) == 0)
    {
      throw std::invalid_argument("window size not supported");
    }
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
  DeriveQuadFirst(int windowSize) : mWindowsize(windowSize)
  {
    if (convolutionWeights.count(windowSize) == 0)
    {
      throw std::invalid_argument("window size not supported");
    }
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

class SmoothGaussian
{
public:
  SmoothGaussian(int windowSize) : mWindowsize(windowSize)
  {
    if (convolutionWeights.count(windowSize) == 0)
    {
      throw std::invalid_argument("window size not supported");
    }
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

      {{{5, {16, 1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0}},
        {7, {64, 1, 6, 15, 20, 15, 6, 1, 0, 0, 0, 0, 0, 0}},
        {9, {256, 1, 8, 28, 56, 70, 56, 28, 8, 1, 0, 0, 0, 0}},
        {11, {2048, 1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1, 0, 0, 0, 0}},
        {13, {8192, 1, 12, 66, 220, 495, 792, 924, 792, 495, 220, 66, 12, 1, 0, 0}},
        {15, {32768, 1, 14, 91, 364, 1001, 2002, 3003, 3432, 3003, 2002, 1001, 364, 91, 14, 1}},
        {17, {131072, 1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820, 560, 120, 16, 1}},
        {19, {524288, 1, 18, 153, 816, 3060, 8568, 18564, 31824, 43758, 48620, 43758, 31824, 18564, 8568, 3060, 816, 153, 18, 1}},
        {21, {2097152, 1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756, 167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1}},
        {23, {8388608, 1, 22, 231, 1540, 7315, 26334, 74613, 170544, 319770, 497420, 646646, 705432, 646646, 497420, 319770, 170544, 74613, 26334, 7315, 1540, 231, 22, 1}},
        {25, {33554432, 1, 24, 276, 2024, 10626, 42504, 134596, 346104, 735471, 1307504, 1961256, 2496144, 2704156, 2496144, 1961256, 1307504, 735471, 346104, 134596, 42504, 10626, 2024, 276, 24, 1}}}}

  };
};

class SmoothAverage
{
public:
  SmoothAverage(int windowSize) : mWindowsize(windowSize)
  {
    if (convolutionWeights.count(windowSize) == 0)
    {
      throw std::invalid_argument("window size not supported");
    }
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

      {{{5, {5, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}},
        {7, {7, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}},
        {9, {9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0}},
        {11, {11, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0}},
        {13, {13, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0}},
        {15, {15, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
        {17, {17, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
        {19, {19, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
        {21, {21, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
        {23, {23, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
        {25, {25, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}}}

  };
};

} // namespace filter
