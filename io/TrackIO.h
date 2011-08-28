#ifndef TRACKIO_H_
#define TRACKIO_H_
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include <hdf5.h>

#include "base/StringUtil.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "base/Common.h"
#include "base/RefCount.h"
#include "base/Hdf5Util.h"
#include "common/Track.h"
#include "io/Traits.h"
#include "common/TrackMetadata.h"
#include "base/Types.h"

class TrackFile;


// A trackfile can store multiple subtracks
// For example, in genomic data these subtracks would correspond to individual chromosomes
class TrackFile
{
public:
    
    typedef boost::shared_ptr<TrackFile> Ptr;
    
    TrackFile()
    : h5file_(-1)
    {}
    
    TrackFile(const std::string& name) 
    : h5file_(-1)
    { Open(name); }
    
    virtual ~TrackFile() {
        Close();
    }
    
    
    // General Track stuff
    DataTypeEnum GetTrackType(const std::string& trackname);
    std::vector<std::string> GetSubTrackNames(const std::string& subtrackname) const;
    std::vector<std::string> GetTrackNames() const;
    TrackMetadata GetSubTrackMetadata(const std::string& trackname, const std::string& subtrackname) const;
    
    bool HasSubTrack(const std::string& trackname, const std::string& subtrackname) const;
    bool HasTrack(const std::string& trackname) const;
    
    // Write Stuff
    template <typename T>
    bool WriteSubTrack(Track<T>& data);
    
    // Read stuff
    template <typename T>
    bool ReadSubTrack(const std::string& trackname,
                      const std::string& subtrackname,
                      Track<T>& subtrack);
    
    bool Open(const char* fname);
    bool Open(const std::string& fname);
    void Close();
    bool IsOpen() { return isopen_; }
    
private:
    DISALLOW_COPY_AND_ASSIGN(TrackFile)
    
    bool Create(const char* fname);
    bool Create(const std::string& fname);
    
    hid_t h5file_;
    std::string filename_;
    bool isopen_;
};

// Implementation
// FUCK YOU C++

//-------------------------------------------------------------------
// Implementation

template<typename DataT>
bool TrackFile::WriteSubTrack(Track<DataT>& subtrack)
{
    hid_t dataset;
    hid_t track_group;
    hid_t root_group;
    hid_t space;
    //  hid_t status;
    std::vector<std::string> tracknames;
    hid_t dcpl;
    
    hsize_t best_chunk_size = 4096*16;
    hsize_t memsize[1] = {(hsize_t) subtrack.size()};
    //const hsize_t chunk_size = std::min(best_chunk_size, memsize[0] / 2);
    //hsize_t chunk[1] = {chunk_size};
    //DataT buff[chunk_size];
    
    if (!isopen_) {
        ERRORLOG("File not open " + filename_);
        return false;    
    }
    
    
    root_group = H5Gopen2(h5file_, "/", H5P_DEFAULT);
    
    // Make sure track exists.
    tracknames = GetTrackNames();
    if (std::find(tracknames.begin(), tracknames.end(), subtrack.trackname()) == tracknames.end()) {
        track_group = H5Gcreate2(root_group, subtrack.trackname().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    } else {
        track_group = H5Gopen2(root_group, subtrack.trackname().c_str(), H5P_DEFAULT);
    }
    if (track_group < 0) {
        ERRORLOG("Can't open track group " + subtrack.trackname());
        H5Gclose(root_group);
        return false;
    }
    
    // Create dataspace, setting max size to NULL so max size = curr size
    space = H5Screate_simple(1, memsize, NULL);
    if (space < 0) {
        ERRORLOG("Can't create data space");
        return false;
    }      
    
    //H5Sset_extent_simple(data_space.id(), 1, memsize, NULL);
    // Check if subtrack exists already
    tracknames = GetSubTrackNames(subtrack.trackname());
    if (std::find(tracknames.begin(), tracknames.end(), subtrack.subtrackname()) == tracknames.end()) {
        
        // Create dataset creation plist, and set chunk size
        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        //status = H5Pset_chunk(dcpl, 1, chunk);
        if (DataTypeTraits<DataT>::IsCompound()) {
            hid_t type = DataTypeTraits<DataT>::H5Type();
            assert(type > 0);
            dataset = H5Dcreate2(track_group, subtrack.subtrackname().c_str(),
                                 type, space,
                                 H5P_DEFAULT, dcpl, H5P_DEFAULT);
            H5Tclose(type);
        } else {
            dataset = H5Dcreate2(track_group, subtrack.subtrackname().c_str(),
                                 DataTypeTraits<DataT>::H5Type(), space, H5P_DEFAULT,
                                 dcpl, H5P_DEFAULT);
        }
    } else {
        
        dataset = H5Dopen2(track_group, subtrack.subtrackname().c_str(), H5P_DEFAULT);
    }
    
    if (dataset < 0) {
        ERRORLOG("Couldn't open subtrack " + subtrack.subtrackname());
        H5Gclose(root_group);
        H5Gclose(track_group);
        return false;
    }
    
    // Write data in chunks, using buff
    /*
     int num_chunks = ceil(((float) memsize[0]) / chunk_size);
     int pos;
     hid_t offset[1];
     hid_t count[1];
     
     
     for (int i = 0; i < num_chunks-1; ++i) {
     for (int j = 0; j < chunk_size; ++j) {
     buff[j] = track.get(pos);
     pos++;
     }
     offset[0] = i;
     status[0] = chunk_size;
     status = H5Sselect_hyperslab(dataset, H5S_SELECT_SET, offset, NULL, count, NULL);
     if (status < 0) {
     ERRORLOG("Couldn't select slab " + i);
     H5Gclose(root_group);
     H5Dclose(dataset);
     H5Gclose(track_group);
     return false;
     }
     status = H5Dwrite(dataset, DataTypeTraits<DataT>::H5Type(),
     space, 
     }*/
    hid_t err = H5Dwrite(dataset, DataTypeTraits<DataT>::H5Type(), H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*subtrack.begin()));
    
    if (err < 0) {
        ERRORLOG("Error writing subtrack");
        H5Gclose(root_group);
        H5Dclose(dataset);
        H5Gclose(track_group);
        return false;
    }
    
    bool attrwritesuccess = true;
    // Write size attributes
    if (!WriteAttribute(dataset, "Start", 1, subtrack.start())) {
        attrwritesuccess = false;
    }
    
    if (!WriteAttribute(dataset, "Stop", 1, subtrack.stop())) {
        attrwritesuccess = false;
    }
    
    if (!WriteAttribute(dataset, "Name", subtrack.subtrackname())) {
        attrwritesuccess = false;
    }
    
    if (!WriteAttribute(dataset, "Resolution", 1, subtrack.resolution())) {
        attrwritesuccess = false;
    }
    if (!WriteAttribute(dataset, "AStart", 1, subtrack.astart())) {
        attrwritesuccess = false;
    }
    if (!WriteAttribute(dataset, "AStop", 1, subtrack.astop())) {
        attrwritesuccess = false;
    }
    if (!attrwritesuccess) {
        ERRORLOG("Couldn't write attributes");
        H5Gclose(root_group);
        H5Gclose(track_group);
        H5Dclose(dataset);
        return false;
    }
    
    H5Gclose(root_group);
    H5Dclose(dataset);
    H5Gclose(track_group);
    return true;
}

template <typename DataT>
bool
TrackFile::ReadSubTrack(const std::string& trackname,
                        const std::string& subtrackname,
                        Track<DataT>& subtrack)
{
    std::string fname = filename_;
    hid_t dataset;
    int start;
    int stop;
    int astart;
    int astop;
    int resolution;
    
    if(!isopen_) {
        ERRORLOG(fname + " is not open");
        return false;
    }
    ScopedH5GOpen root_group(h5file_, std::string("/"));
    // Do some error checking...
    ScopedH5GOpen track_group(root_group.id(), trackname);
    //= H5Gopen2(root_group, trackname.c_str(), H5P_DEFAULT);
    if (track_group < 0) {
        ERRORLOG("Can't open track group " + trackname);
        return false;
    }
    // CHeck if dataset exists before trying to open
    if (!HasSubTrack(trackname, subtrackname)) {
        WARNLOG("Can't find subtrack " + subtrackname);
        return false;
    }
    dataset = H5Dopen2(track_group, subtrackname.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
        ERRORLOG("Can't open subtrack " + subtrackname);
        return false;
    }
    
    // Read size attributes...
    if (!ReadAttribute(dataset, "Start", 1, &start)
        || !ReadAttribute(dataset, "Stop", 1, &stop)
        || !ReadAttribute(dataset, "Resolution", 1, &resolution)
        || !ReadAttribute(dataset, "AStart", 1, &astart)
        || !ReadAttribute(dataset, "AStop", 1, &astop)) {
        ERRORLOG("Can't read attributes");
        H5Dclose(dataset);
        return false;
    }
    subtrack.set_abs_extends(astart, astop);
    subtrack.set_resolution(resolution);
    subtrack.set_extends(start, stop);
    subtrack.set_trackname(trackname);
    subtrack.set_subtrackname(subtrackname);
    herr_t status;
    
    hid_t type = DataTypeTraits<DataT>::H5Type();
    if (DataTypeTraits<DataT>::IsCompound()) {
        status = H5Dread(dataset, type, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, &(*subtrack.begin()));
        H5Tclose(type);
    }  else {
        status  = H5Dread(dataset, DataTypeTraits<DataT>::H5Type(), H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, &(*subtrack.begin()));
    }
    if (status < 0) {
        ERRORLOG("Can't read dataset " + subtrackname);
        H5Dclose(dataset);
        return false;  
    }
    H5Dclose(dataset);
    return true;
}



#endif
