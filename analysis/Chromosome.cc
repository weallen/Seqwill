#include "analysis/Chromosome.h"

//
// CHROMOSOME
//
void Chromosome::Open() { 
  string fname = m_dirname + "/" + m_chrname + ".h5";
  FILE* f = fopen(fname.c_str(), "r");

  if (f == NULL) {
    m_h5file = H5Fcreate(fname.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    assert(m_h5file >= 0);
    H5Gcreate(m_h5file, "/data", 0);
    H5Fclose(m_h5file);
  } else {
    fclose(f);
  }
  m_h5file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(m_h5file >= 0);
  m_isopen = true;
}

int Chromosome::Close() { 
  int ret = H5Fclose(m_h5file); 
  m_isopen = false;
  return ret;
}

int Chromosome::WriteSeq(const string& seq) 
{
  hid_t root, attrspace, startAttr, endAttr, dset, dspace;
  hsize_t dims; 
  int start, end;

  root = H5Gopen(m_h5file, "/");
  attrspace = H5Screate(H5S_SCALAR);
  start = 0;
  end = seq.size();

  startAttr = H5Acreate(root, "start", H5T_NATIVE_INT, attrspace, H5P_DEFAULT);
  endAttr = H5Acreate(root, "end", H5T_NATIVE_INT, attrspace, H5P_DEFAULT);
  H5Awrite(startAttr, H5T_NATIVE_INT, &start);
  H5Awrite(startAttr, H5T_NATIVE_INT, &end);

  dims= seq.size();
  dspace = H5Screate_simple(1, &dims, NULL);
  
  dset = H5Dcreate(m_h5file, "/seq", H5T_NATIVE_UINT8, dspace, H5P_DEFAULT);
  int ret = H5Dwrite(dset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, seq.c_str());
 
  H5Aclose(startAttr);
  H5Aclose(endAttr);
  H5Dclose(dset);
  H5Sclose(dspace);
  H5Sclose(attrspace);
  return ret;
}

int Chromosome::ReadSeq(string* seq) 
{
  hid_t dset;
  int len;
  char* data;
  dset = H5Dopen(m_h5file, "/seq");
  len = GetLength();
  seq->resize(len);
  data = new char[len];
  H5Dread(dset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data); 
  for (int i=0; i < len; ++i) {
    (*seq)[i] = data[i];
  }
  delete data;
  H5Dclose(dset);
  return 1;
}

int Chromosome::GetLength() 
{
  hid_t attr, root;
  if (m_len == -1) {
    root = H5Gopen(m_h5file, "/");
    attr = H5Aopen_name(root, "end");
    H5Aread(attr, H5T_NATIVE_INT, &m_len);
    H5Aclose(attr);
  }
  return m_len;
}

// TODO Finish get track names
svec<string> Chromosome::GetTrackNames() {
  hid_t attr, dspace, dtype, root;
  hsize_t dsize, cell_size, num_cells;
  char* attr_data;
  if (m_tracknames.size() == 0) {
    root = H5Gopen(m_h5file, "/");
    assert(root >= 0);
    attr = H5Aopen_name(root, "tracknames");
    assert(attr >= 0);
    dspace = H5Aget_space(attr);
    assert(dspace >= 0);
    assert(H5Sget_simple_extent_dims(dspace, &num_cells, NULL) == 1);
    assert(H5Sclose(dspace) >= 0);
  }
}

int Chromosome::WriteTrack(const string& name, const ublas::vector<float>& v) 
{
  hid_t grp, dspace, dset;

  float arr[v.size()];
  for (int i=0; i < v.size(); i++) {
    arr[i] = v(i);
  }
  
  grp = H5Gopen(m_h5file, "/data");
  hsize_t dims[] = {v.size()};
  dspace = H5Screate_simple(1, dims, NULL);
  dset = H5Dcreate(grp, name.c_str(), H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT);
  int ret = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr);
  H5Dclose(dset);
  H5Gclose(grp);
  H5Sclose(dspace);
  return ret;
}

ublas::vector<float> Chromosome::ReadTrack(const string& name) 
{
  hid_t dset, dspace; 
  hsize_t curr_dims[1];
  hsize_t max_dims[1];
  string dsetname = "/data/" + name;
  dset = H5Dopen(m_h5file, dsetname.c_str());
  dspace = H5Dget_space(dset);
  H5Sget_simple_extent_dims(dspace, curr_dims, max_dims);
  ublas::vector<float> v(max_dims[0]);
  float data[max_dims[0]];
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  for (int i=0; i < v.size(); i++) {
    v(i) = data[i];
  }
  H5Sclose(dspace);
  H5Dclose(dset);
  return v;
}

void Chromosome::DeleteTrack(const string& name) 
{
  
}
