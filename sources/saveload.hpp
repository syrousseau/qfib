/// \file
// --------------------------------------------------------------------------
// This file is part of the reference implementation for the paper
//    QFib: Fast and Efficient Brain Tractogram Compression
//    C. Mercier*, S. Rousseau*, P. Gori, I. Bloch and T. Boubekeur
//    NeuroInformatics 2020
//    DOI: 10.1007/s12021-020-09452-0
//
// All rights reserved. Use of this source code is governed by a 
// MIT license that can be found in the LICENSE file.
// --------------------------------------------------------------------------
#pragma once

#include "compression.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include "bitvalue.hpp"

using namespace std;

//Version number of this version of the software
static const uint8_t versionNumber = 1;

bool isnan(const Vector3f v)
{
    return ((isnan(v.x()) || isnan(v.y())) || isnan(v.z()));
}

bool isinf(const Vector3f v)
{
    return ((isinf(v.x()) || isinf(v.y())) || isinf(v.z()));
}

namespace fls
{

///
/// \brief loadTckFibers Function to load fibers from a tck file
/// \param filename Path of the file to load
/// \param fibers Table in which the loaded bundle is stored
///

void loadTckFibers(const string & filename, std::vector<std::vector<Vector3f> > & fibers)
{
    fibers.clear();
    if (filename.find(".tck")==string::npos)
    {
        cerr << "File must be in the tck format" << endl;
        exit(EXIT_FAILURE);
    }
    fstream file;
    file.open(filename, ios_base::in | ios_base::binary);
    if (!file.is_open())
    {
        cerr << "Impossible to find file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    string message;
    file >> message;
    if (message!="mrtrix")
    {
        cerr << "Tck file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    file >> message;
    if (message!="tracks")
    {
        cerr << "Tck file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    unsigned offset=0;
    bool count = false, float32le = false;
    while (message!="END")
    {
        file >> message;
        if (message=="count:")
        {
            file >> message;
            cout << stoi(message) << " fibers found ..." << endl;
            fibers.reserve(stoi(message));
            count = true;
        }
        if (message.find("datatype")!=string::npos)
        {
            file >> message;
            if (message!="Float32LE")
            {
                cerr << "Format not supported: " << message << endl;
                exit(EXIT_FAILURE);
            }
            else
                float32le = true;
        }
        if (message.find("file")!=string::npos)
        {
            file >> message;//"."
            file >> message;
            offset = stoi(message);
        }
        if (message.find("quantification")!=string::npos)
        {
            file >> message;
            if ((message=="octahedral" && quantizationMethod != QuantizationMethod::OCTAHEDRAL) ||
                    (message=="fibonacci" && quantizationMethod != QuantizationMethod::SPHERICALFIBONACCI))
            {
                cerr << endl << "----------!WARNING!----------" << endl;
                cerr << "Previous quantification method (" << message << ") was not the one selected here" << endl;
                cerr << "----------!WARNING!----------" << endl << endl;
            }
        }
    }
    if (offset==0)
    {
        cerr << "file needs to be specified in tck file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
    if (!(float32le && count))
    {
        cerr << "Error with file header, count or datatype not specified" << endl;
        exit(EXIT_FAILURE);
    }
    file.seekg(offset, file.beg);
    Vector3f newPoint;
    newPoint << 0,0,0;
    std::vector<Vector3f> newFiber;
    newFiber.clear();
    while (!isinf(newPoint))
    {
        file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
        if (isnan(newPoint) || isinf(newPoint))
        {
            fibers.push_back(newFiber);
            newFiber.clear();
            if (isinf(newPoint))
            {
                break;
            }
            file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
            if (!isnan(newPoint) && !isinf(newPoint))
                newFiber.push_back(newPoint);
        }
        else
        {
            newFiber.push_back(newPoint);
        }
    }
    file.close();
}

///
/// \brief loadCompressAndSave Function to compress fibers out-of-core
/// \param filename_in Path of the input file (tck)
/// \param filename_out Path of the output file (qfib)
///

template<typename T>
void loadCompressAndSave(const string & filename_in, const string & filename_out)
{
    if (filename_in.find(".tck")==string::npos)
    {
        cerr << "File must be in the tck format" << endl;
        exit(EXIT_FAILURE);
    }
    fstream file;
    file.open(filename_in, ios_base::in | ios_base::binary);

    if (filename_out.find(".qfib")==string::npos)
    {
        cerr << "File must be in the qfib format" << endl;
        exit(EXIT_FAILURE);
    }
    fstream file_out;
    file_out.open(filename_out, ios_base::out | ios_base::binary);

    if (!file_out.is_open())
    {
        cerr << "Impossible to find file " << filename_out << endl;
        exit(EXIT_FAILURE);
    }
    string message;
    file >> message;
    if (message!="mrtrix")
    {
        cerr << "Tck file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    file >> message;
    if (message!="tracks")
    {
        cerr << "Tck file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    unsigned nbFibers;
    unsigned offset=0;
    bool count = false, float32le = false;
    while (message!="END")
    {
        file >> message;
        if (message=="count:")
        {
            file >> message;
            cout << "Number of fibers found: " << stoi(message) << endl;
            nbFibers = stoi(message);
            count = true;
        }
        if (message.find("datatype")!=string::npos)
        {
            file >> message;
            if (message!="Float32LE")
            {
                cerr << "Format not supported: " << message << endl;
                exit(EXIT_FAILURE);
            }
            else
                float32le = true;
        }
        if (message.find("file")!=string::npos)
        {
            file >> message;//"."
            file >> message;
            offset = stoi(message);
        }
        if (message.find("quantification")!=string::npos)
        {
            file >> message;
            if ((message=="octahedral" && quantizationMethod != QuantizationMethod::OCTAHEDRAL) ||
                    (message=="fibonacci" && quantizationMethod != QuantizationMethod::SPHERICALFIBONACCI))
            {
                cerr << "----------!WARNING!----------" << endl;
                cerr << "Previous quantification method (" << message << ") was not the one selected here" << endl;
                cerr << "----------!WARNING!----------" << endl;
            }
        }
    }
    if (offset==0)
    {
        cerr << "file needs to be specified in tck file: " << filename_in << endl;
        exit(EXIT_FAILURE);
    }
    file.seekg(offset, file.beg);
    if (!(float32le && count))
    {
        cerr << "Error with file header, count or datatype not specified" << endl;
        exit(EXIT_FAILURE);
    }

    float mdot = 0.f; //compute ratio here
    Vector3f newPoint;
    std::vector<Vector3f> newFiber;
    while (!isinf(newPoint))
    {
        file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
        if (isnan(newPoint) || isinf(newPoint))
        {
            mdot = max(fc::computeRatio(newFiber), mdot);
            newFiber.clear();
            if (isinf(newPoint))
            {
                break;
            }
            file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
            if (!isnan(newPoint) && !isinf(newPoint))
                newFiber.push_back(newPoint);
        }
        else
        {
            newFiber.push_back(newPoint);
        }
    }
    epsMapping = 0.f;
    for(unsigned i = 0; i < 20; ++i)
        epsMapping = 3 * sqrt(2) * pow(2, -bitcount/2.f) * ((1 - (mdot))/2+epsMapping);
    //Error maximization proportionaly to itself (epsMapping < 1)
    epsMapping = sqrt(epsMapping);
    float ratio = std::min(1.f, (1.f - mdot) / 2.f + epsMapping);

    file.seekg(offset, file.beg);

    file_out.write(reinterpret_cast<const char *>(&versionNumber), sizeof (uint8_t));
    file_out.write(reinterpret_cast<const char *>(&nbFibers), sizeof (unsigned));
    file_out.write(reinterpret_cast<const char *>(&ratio), sizeof (float));
    uint8_t method;
    if (quantizationMethod == QuantizationMethod::OCTAHEDRAL)
        method = 0;
    else//fibonacci
        method = 1;
    file_out.write(reinterpret_cast<const char *>(&method), sizeof (uint8_t));
    file_out.write(reinterpret_cast<const char *>(&bitcount), sizeof (uint8_t));

    unsigned nbOfThreads = 10000;//omp_get_max_threads();
    unsigned nbThreadsUsed=0;

    newFiber.clear();
    newPoint << 0, 0, 0;
    CompressedFiber<T> newCFiber;
    std::vector<CompressedFiber<T>> fibers(nbOfThreads);
#pragma omp parallel
#pragma omp single
    while (!isinf(newPoint))
    {
        file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
        if (isnan(newPoint) || isinf(newPoint))
        {
#pragma omp task firstprivate(newFiber, newCFiber, nbThreadsUsed)
            {
                fc::compressFiber(newFiber, ratio, newCFiber);
                fibers[nbThreadsUsed]=(newCFiber);
            }
            nbThreadsUsed++;
            if (nbThreadsUsed>=nbOfThreads)
            {
#pragma omp taskwait
                {
                    nbThreadsUsed = 0;
                    for (unsigned i=0; i<nbOfThreads; i++)
                    {
                        Vector3f origin = fibers[i].origin;
                        file_out.write(reinterpret_cast<const char *>(&origin), sizeof (Vector3f));
                        Vector3f second = fibers[i].second;
                        file_out.write(reinterpret_cast<const char *>(&second), sizeof (Vector3f));
                        uint16_t nbPoints = fibers[i].data.size();
                        file_out.write(reinterpret_cast<const char *>(&nbPoints), sizeof (uint16_t));
                        for (unsigned j=0; j<nbPoints; j++)
                        {
                            T value = fibers[i].data[j].getCombinedValue();//qv;
                            file_out.write(reinterpret_cast<const char *>(&value), sizeof (T));
                        }
                        fibers[i].data.clear();
                    }
                }
            }
            newFiber.clear();
            if (isinf(newPoint))
            {
                break;
            }
            file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
            if (!isnan(newPoint) && !isinf(newPoint))
                newFiber.push_back(newPoint);
        }
        else
        {
            newFiber.push_back(newPoint);
        }
    }
    for (unsigned i=0; i<nbThreadsUsed; i++)
    {
        Vector3f origin = fibers[i].origin;
        file_out.write(reinterpret_cast<const char *>(&origin), sizeof (Vector3f));
        Vector3f second = fibers[i].second;
        file_out.write(reinterpret_cast<const char *>(&second), sizeof (Vector3f));
        uint16_t nbPoints = fibers[i].data.size();
        file_out.write(reinterpret_cast<const char *>(&nbPoints), sizeof (uint16_t));
        for (unsigned j=0; j<nbPoints; j++)
        {
            T value = fibers[i].data[j].getCombinedValue();
            file.write(reinterpret_cast<const char *>(&value), sizeof (T));
        }
    }
    file.close();
    file_out.close();
}

///
/// \brief loadDecompressAndSave Function to decompress fibers out_of_core
/// \param filename_in Path of the input file (qfib)
/// \param filename_out Path of the output file (tck)
///

template<typename T>
void loadDecompressAndSave(const string & filename_in, const string & filename_out)
{
    fstream file_out;
    if (filename_out.find(".tck")==string::npos)
    {
        cerr << "Error, file needs to be in tck format" << endl;
        exit(EXIT_FAILURE);
    }
    file_out.open(filename_out, ios_base::out | ios_base::binary);
    if (!file_out.is_open())
    {
        cerr << "Error, impossible to write file " << filename_out << endl;
        exit(EXIT_FAILURE);
    }
    fstream file_in;
    if (filename_in.find(".qfib")==string::npos)
    {
        cerr << "Error, file needs to be in qfib format" << endl;
        exit(EXIT_FAILURE);
    }
    file_in.open(filename_in, ios_base::in | ios_base::binary);
    if (!file_in.is_open())
    {
        cerr << "Impossible to read file " << filename_in << endl;
        exit(EXIT_FAILURE);
    }
    uint8_t version;
    file_in.read(reinterpret_cast<char *>(&version), sizeof (uint8_t));
    if (version != versionNumber)
    {
        cerr << endl << "!!! Warning !!" << endl;
        cerr << "Version numbers of file and storage are not the same" << endl;
        cerr << "!!! Warning !!" << endl << endl;
    }
    unsigned nbFibers;
    file_in.read(reinterpret_cast<char *>(&nbFibers), sizeof (unsigned));
    float ratio;
    file_in.read(reinterpret_cast<char *>(&ratio), sizeof (float));
    uint8_t method;
    file_in.read(reinterpret_cast<char *>(&method), sizeof (uint8_t));
    if (method == 0)//octahedral
        quantizationMethod = QuantizationMethod::OCTAHEDRAL;
    else//fibonacci
        quantizationMethod = QuantizationMethod::SPHERICALFIBONACCI;
    uint8_t bits;
    file_in.read(reinterpret_cast<char *>(&bits), sizeof (uint8_t));

    file_out << "mrtrix tracks" << endl;
    file_out << "count: " << nbFibers << endl;
    file_out << "datatype: Float32LE" << endl;
    file_out << "quantification: ";
    if (quantizationMethod == QuantizationMethod::OCTAHEDRAL)
        file_out << "octahedral" << endl;
    else
        file_out << "fibonacci" << endl;
    uint64_t offset = static_cast<unsigned>(file_out.tellp())+15;
    file_out << "file: . " << offset << endl;
    file_out << "END" << endl;
    file_out.seekp(offset, file_out.beg);
    Vector3f nanVec;
    nanVec << nanf(""), nanf(""), nanf("");
    Vector3f infVec;
    infVec << numeric_limits<float>::infinity(), numeric_limits<float>::infinity(), numeric_limits<float>::infinity();

    unsigned nbOfThreads = 10000;
    unsigned nbThreadsUsed=0;
    std::vector<std::vector<Vector3f> > uncompressedFiber(nbOfThreads);

    offset = static_cast<unsigned>(file_out.tellp());
    file_out.close();
#pragma omp parallel
#pragma omp single
    for (unsigned i=0; i<nbFibers; i++)
    {
        CompressedFiber<T> fiberC;
        Vector3f origin;
        file_in.read(reinterpret_cast<char *>(&origin), sizeof (Vector3f));
        fiberC.origin = origin;
        Vector3f second;
        file_in.read(reinterpret_cast<char *>(&second), sizeof (Vector3f));
        fiberC.second = second;
        uint16_t nbPoints;
        file_in.read(reinterpret_cast<char *>(&nbPoints), sizeof (uint16_t));
        fiberC.data.reserve(nbPoints);
        for (unsigned j=0; j<nbPoints; j++)
        {
            CompressedVector<T> cv;
            T value;
            file_in.read(reinterpret_cast<char *>(&value), sizeof (T));
            cv.setCombinedValue(value);
            fiberC.data.push_back(cv);
        }
#pragma omp task firstprivate(fiberC, nbThreadsUsed)
        {
            fc::decompressFiber(fiberC, ratio, uncompressedFiber[nbThreadsUsed]);
        }
        nbThreadsUsed++;
        if (nbThreadsUsed>=nbOfThreads)
        {
#pragma omp taskwait
            {
                //Write fiber in tck file
                nbThreadsUsed=0;
                file_out.open(filename_out);
                file_out.seekp(offset, file_out.beg);
                for (unsigned k=0; k<nbOfThreads; k++)
                {
                    for (unsigned j=0; j<uncompressedFiber[k].size(); j++)
                    {
                        file_out.write(reinterpret_cast<const char *>(&uncompressedFiber[k][j]), sizeof (Vector3f));
                    }
                    uncompressedFiber[k].clear();
                    file_out.write(reinterpret_cast<const char *>(&nanVec), sizeof (Vector3f));
                }
                offset = file_out.tellp();
                file_out.close();
            }
        }
    }

    file_out.open(filename_out);
    file_out.seekp(offset, file_out.beg);
    for (unsigned k=0; k<nbThreadsUsed; k++)
    {
        for (unsigned j=0; j<uncompressedFiber[k].size(); j++)
        {
            file_out.write(reinterpret_cast<const char *>(&uncompressedFiber[k][j]), sizeof (Vector3f));
        }
        uncompressedFiber[k].clear();
        file_out.write(reinterpret_cast<const char *>(&nanVec), sizeof (Vector3f));
    }
    file_out.write(reinterpret_cast<const char *>(&infVec), sizeof (Vector3f));
    file_in.close();
    file_out.close();
}

void saveTck(string filename, const std::vector<std::vector<Vector3f> > & fibers)
{
    fstream file;
    if (filename.find(".tck")==string::npos)
        filename = filename+".tck";
    file.open(filename, ios_base::out | ios_base::binary);
    if (!file.is_open())
    {
        cerr << "Impossible to write file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    file << "mrtrix tracks" << endl;
    file << "count: " << fibers.size() << endl;
    file << "datatype: Float32LE" << endl;
    file << "quantification: ";
    if (quantizationMethod == QuantizationMethod::OCTAHEDRAL)
        file << "octahedral" << endl;
    else
        file << "fibonacci" << endl;
    unsigned offset = static_cast<unsigned>(file.tellp())+15;
    file << "file: . " << offset << endl;
    file << "END" << endl;
    file.seekp(offset, file.beg);
    Vector3f nanVec;
    nanVec << nanf(""), nanf(""), nanf("");
    Vector3f infVec;
    infVec << numeric_limits<float>::infinity(), numeric_limits<float>::infinity(), numeric_limits<float>::infinity();
    for (unsigned i=0; i<fibers.size(); i++)
    {
        for (unsigned j=0; j<fibers[i].size(); j++)
        {
            file.write(reinterpret_cast<const char *>(&fibers[i][j]), sizeof (Vector3f));
        }
        file.write(reinterpret_cast<const char *>(&nanVec), sizeof (Vector3f));
    }
    file.write(reinterpret_cast<const char *>(&infVec), sizeof (Vector3f));
    file.close();
}

///
/// \brief saveCompressedFibers Function to save fibers in the qfib format
/// \param filename Path of the output file (qfib)
/// \param fibers Fibers compressed to save in the file
/// \param ratio Ratio of the compressed bundle (necessary to decompress)
///

template<typename T>
void saveCompressedFibers(const string filename, const std::vector<CompressedFiber<T>> & fibers, const float ratio)
{
    fstream file;
    if (filename.find("qfib")==string::npos)
    {
        cerr << "File to record must be a qfib file" << endl;
        exit(EXIT_FAILURE);
    }
    file.open(filename, ios_base::out | ios_base::binary);
    if (!file.is_open())
    {
        cerr << "Impossible to create file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    file.write(reinterpret_cast<const char *>(&versionNumber), sizeof (uint8_t));
    unsigned nbFibers = fibers.size();
    file.write(reinterpret_cast<const char *>(&nbFibers), sizeof (unsigned));
    file.write(reinterpret_cast<const char *>(&ratio), sizeof (float));
    uint8_t method;
    if (quantizationMethod == QuantizationMethod::OCTAHEDRAL)
        method = 0;
    else//fibonacci
        method = 1;
    file.write(reinterpret_cast<const char *>(&method), sizeof (uint8_t));
    uint8_t bits = sizeof (T) * CHAR_BIT;
    file.write(reinterpret_cast<const char *>(&bits), sizeof (uint8_t));
    for (unsigned i=0; i<nbFibers; i++)
    {
        Vector3f origin = fibers[i].origin;
        file.write(reinterpret_cast<const char *>(&origin), sizeof (Vector3f));
        Vector3f second = fibers[i].second;
        file.write(reinterpret_cast<const char *>(&second), sizeof (Vector3f));
        uint16_t nbPoints = fibers[i].data.size();
        file.write(reinterpret_cast<const char *>(&nbPoints), sizeof (uint16_t));
        for (unsigned j=0; j<nbPoints; j++)
        {
            T value = fibers[i].data[j].getCombinedValue();//qv;
            file.write(reinterpret_cast<const char *>(&value), sizeof (T));
        }
    }
    file.close();
}

///
/// \brief is16bits Function to read from a qfib file whether the file is in 16 bits precision or not (8 bits assumed otherwised)
/// \param filename Path of the qfib file to read
/// \return Return true if the file is saved in 16 bits precision, false otherwise
///

bool is16bits(string filename)
{
    fstream file;
    if (filename.find(".qfib")==string::npos)
    {
        cerr << "File not in qfib format" << endl;
        exit(EXIT_FAILURE);
    }
    file.open(filename, ios_base::in | ios_base::binary);
    if (!file.is_open())
    {
        cerr << "Impossible to read file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    uint8_t version;
    file.read(reinterpret_cast<char *>(&version), sizeof (uint8_t));
    unsigned nbFibers;
    file.read(reinterpret_cast<char *>(&nbFibers), sizeof (unsigned));
    float ratio;
    file.read(reinterpret_cast<char *>(&ratio), sizeof (float));
    uint8_t method;
    file.read(reinterpret_cast<char *>(&method), sizeof (uint8_t));
    if (method==0)
        quantizationMethod = QuantizationMethod::OCTAHEDRAL;
    else
        quantizationMethod = QuantizationMethod::SPHERICALFIBONACCI;
    uint8_t bits;
    file.read(reinterpret_cast<char *>(&bits), sizeof (uint8_t));
    file.close();
    return (bits == 16);
}

///
/// \brief loadCompressedFibers Function to load compressed fibers from a qfib file
/// \param filename Path of the qfib file to read
/// \param fibers Empty vector of compressed fibers in which fibers will be saved
/// \param ratio Ratio read from the file (useful for decompression)
///

template<typename T>
void loadCompressedFibers(string filename, vector<CompressedFiber<T>> & fibers, float& ratio)
{
    fstream file;
    if (filename.find(".qfib")==string::npos)
    {
        cerr << "File not in qfib format" << endl;
        exit(EXIT_FAILURE);
    }
    file.open(filename, ios_base::in | ios_base::binary);
    if (!file.is_open())
    {
        cerr << "Impossible to read file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    uint8_t version;
    file.read(reinterpret_cast<char *>(&version), sizeof (uint8_t));
    if (version != versionNumber)
    {
        cerr << endl << "!!! Warning !!" << endl;
        cerr << "Version numbers of file and storage are not the same" << endl;
        cerr << "!!! Warning !!" << endl << endl;
    }
    unsigned nbFibers;
    file.read(reinterpret_cast<char *>(&nbFibers), sizeof (unsigned));
    fibers.clear();
    fibers.reserve(nbFibers);
    file.read(reinterpret_cast<char *>(&ratio), sizeof (float));
    uint8_t method;
    file.read(reinterpret_cast<char *>(&method), sizeof (uint8_t));
    if (method == 0)//octahedral
        quantizationMethod = QuantizationMethod::OCTAHEDRAL;
    else//fibonacci
        quantizationMethod = QuantizationMethod::SPHERICALFIBONACCI;
    uint8_t bits;
    file.read(reinterpret_cast<char *>(&bits), sizeof (uint8_t));
    for (unsigned i=0; i<nbFibers; i++)
    {
        CompressedFiber<T> fiberC;
        Vector3f origin;
        file.read(reinterpret_cast<char *>(&origin), sizeof (Vector3f));
        fiberC.origin = origin;
        Vector3f second;
        file.read(reinterpret_cast<char *>(&second), sizeof (Vector3f));
        fiberC.second = second;
        uint16_t nbPoints;
        file.read(reinterpret_cast<char *>(&nbPoints), sizeof (uint16_t));
        fiberC.data.reserve(nbPoints);
        for (unsigned j=0; j<nbPoints; j++)
        {
            CompressedVector<T> cv;
            T value;
            file.read(reinterpret_cast<char *>(&value), sizeof (T));
            cv.setCombinedValue(value);
            fiberC.data.push_back(cv);
        }
        fibers.push_back(fiberC);
    }
    file.close();
}

///
/// \brief loadFAValues Function to load FA values from a mif file (MRtrix format)
/// \param filename File in mif format containing the FA values
/// \param valueGrid Vector in which values will be stored (A value in each grid cell)
/// \param dim Empty Vector in which the dimensions of the grid will be stored
/// \param transform Transformation matrix read from file
/// \param vox Vector in which the voxel size in mm is saved
///

void loadFAValues(const string filename, std::vector<float> & valueGrid, std::vector<unsigned> & dim, MatrixXf & transform, std::vector<float> & vox)
{
    valueGrid.clear();
    if (filename.find(".mif")==string::npos)
    {
        cerr << "FA file must be in the mif format" << endl;
        exit(EXIT_FAILURE);
    }
    fstream file;
    file.open(filename, ios_base::in | ios_base::binary);
    if (!file.is_open())
    {
        cerr << "Impossible to find file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    string message;
    file >> message;
    if (message!="mrtrix")
    {
        cerr << "Mif file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    file >> message;
    if (message!="image")
    {
        cerr << "Mif file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    unsigned offset=0;
    bool count = false, float32le = false;
    transform.resize(3,4);
    std::vector<unsigned> layout;
    std::vector<bool> order;
    unsigned position = 0;
    while (message!="END")
    {
        file >> message;
        if (message=="dim:")
        {
            file >> message;
            unsigned pos = message.find_first_of(",");
            unsigned lastpos = message.find_last_of(",");
            dim.push_back(stoi(message.substr(0, pos)));
            dim.push_back(stoi(message.substr(pos+1, lastpos-pos-1)));
            dim.push_back(stoi(message.substr(lastpos+1, message.size())));
        }
        if (message=="vox:")
        {
            file >> message;
            unsigned pos = message.find_first_of(",");
            unsigned lastpos = message.find_last_of(",");
            vox.push_back(stof(message.substr(0, pos)));
            vox.push_back(stof(message.substr(pos+1, lastpos-pos-1)));
            vox.push_back(stof(message.substr(lastpos+1, message.size())));
        }
        if (message=="layout:")
        {
            file >> message;
            unsigned pos = message.find_first_of(",");
            unsigned lastpos = message.find_last_of(",");
            order.push_back(message.front()=='+');
            order.push_back(message.at(pos+1)=='+');
            order.push_back(message.at(lastpos+1)=='+');
            layout.push_back(stoi(message.substr(1, pos)));
            layout.push_back(stoi(message.substr(pos+2, lastpos-pos-2)));
            layout.push_back(stoi(message.substr(lastpos+2, message.size())));
        }
        if (message=="transform:")
        {
            for (unsigned i=0; i<3; i++)
            {
                file >> message;
                transform(position, i) = stof(message.substr(0, message.size()-1));
            }

            file >> message;
            transform(position, 3) = stof(message.substr(0, message.size()));
            position++;
        }
        if (message.find("datatype")!=string::npos)
        {
            file >> message;
            if (message!="Float32LE")
            {
                cerr << "Format not supported: " << message << endl;
                exit(EXIT_FAILURE);
            }
            else
                float32le = true;
        }
        if (message.find("file:")!=string::npos)
        {
            file >> message;//"."
            file >> message;
            offset = stoi(message);
        }
    }
    if (offset==0)
    {
        cerr << "file needs to be specified in mif file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
    file.seekg(offset, file.beg);
    std::vector<unsigned> dimIter(3);
    for (unsigned i=0; i<3; i++)
        dimIter[i] = dim[layout[i]];
    for (unsigned i=0; i<3; i++)
        dim[i] = dimIter[i];
    valueGrid.resize(dimIter[0]*dimIter[1]*dimIter[2]);
    std::vector<unsigned> c(3);
    float value;
    for (unsigned k=0; k<dimIter[2]; k++)
    {
        for (unsigned j=0; j<dimIter[1]; j++)
        {
            for (unsigned i=0; i<dimIter[0]; i++)
            {
                file.read(reinterpret_cast<char *>(&value), sizeof (float));
                if (order[0])
                    c[0] = i;
                else
                    c[0] = dimIter[0]-i;
                if (order[1])
                    c[1] = j;
                else
                    c[1] = dimIter[1]-j;
                if (order[2])
                    c[2] = k;
                else
                    c[2] = dimIter[2]-k;
                valueGrid[c[0]+c[1]*dimIter[0]+c[2]*dimIter[0]*dimIter[1]] = value;
            }
        }
    }
    file.close();
}

///
/// \brief saveDifferenceAsVTK Function to save difference between two tck files in a vtk file containing fibers of the first files and colors depending on the error in mm per point
/// Both input files need to have the exact same number of points and fibers
/// This function is here to compare original fibers to fibers having being compressed and decompressed
/// \param filename1 Path of the input file 1 (tck)
/// \param filename2 Path of the input file 2 (tck)
///

void saveDifferenceAsVTK(string filename1, string filename2)
{
    std::vector<std::vector<Vector3f>> fibers1;
    std::vector<std::vector<Vector3f>> fibers2;
    loadTckFibers(filename1, fibers1);
    loadTckFibers(filename2, fibers2);

    if(fibers1.size() != fibers2.size())
    {
        cout << fibers1.size() << " " << fibers2.size() << endl;
        cerr << "Error, files must contain the same number of fibers" << endl;
        exit(EXIT_FAILURE);
    }
    unsigned nbFibers = fibers1.size();
    unsigned nbPoints = 0;

    for (unsigned i=0; i<nbFibers; i++)
        nbPoints += fibers1[i].size();

    fstream file;
    string filename_out = filename2 + ".vtk";
    file.open(filename_out, ios_base::out);

    file << "# vtk DataFile Version 7.1" << endl;
    file << "Difference file between compressed and decompressed tractograms, obtained with qfib made by Corentin Mercier and Sylvain Rousseau" << endl;
    file << "ASCII" << endl;
    file << "DATASET POLYDATA" << endl;
    file << "POINTS " << nbPoints <<  " float" << endl;

    for (unsigned i=0; i<nbFibers; i++)
    {
        for (unsigned j=0; j<fibers1[i].size(); j++)
            file << fibers1[i][j].x() << " " << fibers1[i][j].y() << " " << fibers1[i][j].z() << endl;
    }

    file << "LINES " << nbFibers << " " << nbPoints+nbFibers << endl;
    uint64_t counter = 0;
    for (unsigned i=0; i<nbFibers; i++)
    {
        file << fibers1[i].size() << " ";
        for (unsigned j=0; j < fibers1[i].size(); j++)
        {
            file << counter << " ";
            counter++;
        }
        file << endl;
    }
    file << endl;
    file << "POINT_DATA " << counter << endl;
    file << "SCALARS acs float 1" << endl;
    file << "LOOKUP_TABLE default" << endl;
    //Computation of the maximum distance
    //        float maxDistance = 0;
    //        for (unsigned i=0; i<fibers1.size(); i++)
    //            for (unsigned j=0; j<fibers1[i].size(); j++)
    //                maxDistance = std::max(maxDistance, (fibers1[i][j] - fibers2[i][j]).norm());

    unsigned worstFiber = 0;
    float ecartMax = 0;
    for (unsigned i=0; i<fibers1.size(); i++)
        for (unsigned j=0; j<fibers1[i].size(); j++)
        {
            float value = (fibers1[i][j] - fibers2[i][j]).norm();
            if (value > ecartMax)
            {
                ecartMax = value;
                worstFiber = i;
            }
            file << value << endl;
        }
    file.close();

    //File with worst fiber : 1st one = original, 2nd one = compressed/decompressed

    string filename_out2 = filename2 + "worst.vtk";
    fstream file2;
    file2.open(filename_out2, ios_base::out);

    file2 << "# vtk DataFile Version 7.1" << endl;
    file2 << "Difference file between compressed and decompressed tractograms, obtained with qfib made by Corentin Mercier and Sylvain Rousseau" << endl;
    file2 << "ASCII" << endl;
    file2 << "DATASET POLYDATA" << endl;
    unsigned f1 = 35163; //Number of the first fiber to draw
    unsigned f2 = 48573; //Number of the second fiber to draw
    //Number of fibers should be written manually and be the worst case of different cases (precision, quantization...)
    unsigned worstNbPoints = fibers1[f1].size() + fibers1[f2].size();
    file2 << "POINTS " << worstNbPoints*2 <<  " float" << endl;

    for (unsigned j=0; j<fibers1[f1].size(); j++)
        file2 << fibers1[f1][j].x() << " " << fibers1[f1][j].y() << " " << fibers1[f1][j].z() << endl;
    for (unsigned j=0; j<fibers1[f1].size(); j++)
        file2 << fibers2[f1][j].x() << " " << fibers2[f1][j].y() << " " << fibers2[f1][j].z() << endl;
    for (unsigned j=0; j<fibers1[f2].size(); j++)
        file2 << fibers1[f2][j].x() << " " << fibers1[f2][j].y() << " " << fibers1[f2][j].z() << endl;
    for (unsigned j=0; j<fibers1[f2].size(); j++)
        file2 << fibers2[f2][j].x() << " " << fibers2[f2][j].y() << " " << fibers2[f2][j].z() << endl;

    file2 << "LINES 4 " << worstNbPoints*2+4 << endl;
    counter = 0;

    //cout << "Worst fiber is number: " << worstFiber << endl;

    for (unsigned i=0; i<2; i++)
    {
    file2 << fibers1[f1].size() << " ";
    for (unsigned j=0; j < fibers1[f1].size(); j++)
    {
        file2 << counter << " ";
        counter++;
    }
    }
    for (unsigned i=0; i<2; i++)
    {
    file2 << fibers1[f2].size() << " ";
    for (unsigned j=0; j < fibers1[f2].size(); j++)
    {
        file2 << counter << " ";
        counter++;
    }
    }
    file2 << endl;

    file2 << endl;
    file2 << "POINT_DATA " << counter << endl;
    file2 << "SCALARS acs float 1" << endl;
    file2 << "LOOKUP_TABLE default" << endl;

    for (unsigned j=0; j < fibers1[f1].size(); j++)
        file2 << "0" << endl;
    for (unsigned j=0; j < fibers2[f1].size(); j++)
        file2 << (fibers1[f1][j] - fibers2[f1][j]).norm() << endl;
    for (unsigned j=0; j < fibers1[f2].size(); j++)
        file2 << "0" << endl;
    for (unsigned j=0; j < fibers2[f2].size(); j++)
        file2 << (fibers1[f2][j] - fibers2[f2][j]).norm() << endl;

    file2.close();
}
}
