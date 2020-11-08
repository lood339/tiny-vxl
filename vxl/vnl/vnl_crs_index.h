// This is core/vnl/vnl_crs_index.h
#ifndef vnl_crs_index_h_
#define vnl_crs_index_h_
//:
// \file
// \brief Compressed Row Storage (CRS) indexing
// \author Matt Leotta (Brown)
// \date   April 13, 2005
//
// \verbatim
//  Modifications
// \endverbatim
//
#include <vector>
#include <utility>
#include <vnl/vnl_export.h>

//: Represents the configuration of a sparse matrix but not the data
//  This is essentially a sparse matrix of indices into a data vector
//  Compressed row storage is used for representation
//  This class is useful when working with several sparse matrices that
//  share a common sparse structure.
class VNL_EXPORT vnl_crs_index
{
 public:
  typedef std::pair<int,int> idx_pair;
  typedef std::vector<idx_pair> sparse_vector;

  //: Constructor - default
  vnl_crs_index() : col_idx_(), row_ptr_() {}

  //: Constructor - from a binary mask
  vnl_crs_index(const std::vector<std::vector<bool> >& mask);

  //: Destructor
  ~vnl_crs_index()= default;

  //: number of rows in the sparse matrix
  int num_rows() const { return int(row_ptr_.size())-1; }

  //: number of columns in the sparse matrix
  int num_cols() const { return num_cols_; }

  //: number of non-zero elements
  int num_non_zero() const { return int(col_idx_.size()); }

  //: returns row \p i as a vector of index-column pairs
  sparse_vector sparse_row(int i) const;

  //: returns column \p j as a vector of index-row pairs
  // \note because of CRS this method is a bit less efficient than sparse_row
  sparse_vector sparse_col(int j) const;

  //: return the index at location (i,j)
  //  returns -1 if the entry is 0
  int operator() (int i, int j) const;

 private:
  //: The number of columns in the matrix
   unsigned int num_cols_{0};
   //: The column for each non-zero element
   std::vector<int> col_idx_;
   //: The index of the first non-zero element in each row
   std::vector<int> row_ptr_;
};

// copy from .cpp
//: Constructor - from a binary mask
vnl_crs_index::vnl_crs_index(const std::vector<std::vector<bool>> & mask)
: num_cols_(mask[0].size())
, col_idx_()
, row_ptr_(mask.size() + 1, 0)
{
    int k = 0;
    for (unsigned int i = 0; i < mask.size(); ++i)
    {
        const std::vector<bool> & col = mask[i];
        row_ptr_[i] = k;
        for (unsigned int j = 0; j < num_cols_; ++j)
        {
            if (col[j])
            {
                col_idx_.push_back(j);
                ++k;
            }
        }
    }
    row_ptr_[mask.size()] = k;
}


//: return the index at location (i,j)
//  returns -1 if the entry is 0
int
vnl_crs_index::operator()(int i, int j) const
{
    int low = row_ptr_[i];
    int high = row_ptr_[i + 1] - 1;
    
    // binary search for finding the element at column j
    while (low <= high)
    {
        if (j < col_idx_[low] || j > col_idx_[high])
            return -1; // element is zero (no index)
        
        int mid = (low + high) >> 1; //(low+high)/2;
        if (j < (int)col_idx_[mid])
            high = mid - 1;
        else if (j > (int)col_idx_[mid])
            low = mid + 1;
        else
            return mid;
    }
    
    return -1; // element is zero (no index)
}


//: returns row \p i as a vector of index-column pairs
vnl_crs_index::sparse_vector
vnl_crs_index::sparse_row(int i) const
{
    sparse_vector row;
    for (int j = row_ptr_[i]; j < row_ptr_[i + 1]; ++j)
    {
        row.push_back(idx_pair(j, col_idx_[j]));
    }
    return row;
}


//: returns column \p j as a vector of index-row pairs
// \note because of CRS this method is a bit less efficient than sparse_row
vnl_crs_index::sparse_vector
vnl_crs_index::sparse_col(int j) const
{
    sparse_vector col;
    for (int i = 0; i < num_rows(); ++i)
    {
        int idx = (*this)(i, j);
        if (idx >= 0)
            col.push_back(idx_pair(idx, i));
    }
    
    return col;
}

#endif // vnl_crs_index_h_
