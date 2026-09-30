R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void insertion_sort(float *arr, const int n)
{
  for (int i = 1; i < n; i++)
  {
    float key = arr[i];
    int   j = i - 1;
    while ((j >= 0) && (key < arr[j]))
    {
      arr[j + 1] = arr[j];
      j--;
    }
    arr[j + 1] = key;
  }
}
)""
