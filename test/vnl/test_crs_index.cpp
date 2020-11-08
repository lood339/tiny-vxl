#include <iostream>
#include <vnl/vnl_crs_index.h>

#include <gtest/gtest.h>

void
display_mask(const std::vector<std::vector<bool>> & mask)
{
  for (const auto & i : mask)
  {
    for (unsigned int j = 0; j < i.size(); ++j)
    {
      std::cout << (i[j] ? '1' : '0') << ' ';
    }
    std::cout << '\n';
  }
  std::cout << std::endl;
}

TEST(vnl_crs_index, test_crs_index)
{
  std::vector<bool> null_col(8, false);
  std::vector<std::vector<bool>> mask(10, null_col);

  mask[0][1] = true;
  mask[0][2] = true;
  mask[0][4] = true;
  mask[2][0] = true;
  mask[2][1] = true;
  mask[2][7] = true;
  mask[3][4] = true;
  mask[4][3] = true;
  mask[4][6] = true;
  mask[6][0] = true;
  mask[6][2] = true;
  mask[6][3] = true;
  mask[6][4] = true;
  mask[6][6] = true;
  mask[6][7] = true;
  mask[7][7] = true;
  mask[9][6] = true;

  display_mask(mask);

  vnl_crs_index crs(mask);

  EXPECT_EQ(crs.num_rows(), (int)mask.size());
  EXPECT_EQ(crs.num_cols(), (int)mask[0].size());
  EXPECT_EQ(crs.num_non_zero(), 17);

  // look up each index and check it's validity
  bool valid = true;
  int curr = 0;
  for (int i = 0; i < crs.num_rows(); ++i)
  {
    for (int j = 0; j < crs.num_cols(); ++j)
    {
      int idx = crs(i, j);
      if (((idx < 0) == mask[i][j]) && idx != curr++)
        valid = false;
      std::cout << crs(i, j) << ' ';
    }
    std::cout << std::endl;
  }
  EXPECT_EQ(valid, true)<<"operator ()\n";

  // test column extraction
  vnl_crs_index::sparse_vector col = crs.sparse_col(1);
  valid = true;
  for (auto & i : col)
  {
    if (i.first != crs(i.second, 1))
      valid = false;
  }
  EXPECT_EQ(col.size() == 2 && valid, true)<<"sparse_col\n";
  col = crs.sparse_col(5);
  EXPECT_EQ(col.size(), 0)<<"empty sparse_col\n";

  // test row extraction
  vnl_crs_index::sparse_vector row = crs.sparse_row(6);
  valid = true;
  for (auto & i : row)
  {
    if (i.first != crs(6, i.second))
      valid = false;
  }
  EXPECT_EQ(row.size() == 6 && valid, true)<<"sparse_row\n";
  row = crs.sparse_row(5);
  EXPECT_EQ(row.size(), 0)<<"empty sparse_row\n";
}


