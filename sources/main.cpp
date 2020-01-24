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
//#define OMP_THREAD_LIMIT 1;

#include <cassert>
#include <iostream>

#include "timer.h"
#include "compression.hpp"
#include "utilities.hpp"
#include "saveload.hpp"
#include "values.hpp"
#include <omp.h>
#include <getopt.h>
#include <sys/stat.h>

///
/// \brief GetFileSize Function to obtain file size
/// \param filename File to consider
/// \return Size of the file
///

long GetFileSize(const std::string & filename)
{
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

///
/// \brief compress Compress a set of fibers
/// \param inputFileName Path of the input file (tck)
/// \param outputFileName Path of the output file (qfib)
/// \param outOfCore Boolean : true when outOfCore, false otherwise
/// \param timings When true, display the compression timings
/// \param displayError When true, compute the compression error and diplay it - only in in-core mode
/// \param displayDetailedError When true, compute the compression detailed error per fiber length and display it - only in in-core mode
/// \param verbose When true, display the maximum angle of the set
/// \param verif When true, verify that segments length are all the same length with a 10% margin - only in in-core mode
/// Deactivate for full performance
/// \param fa Path of an fa file to use for comparison - only in in-core mode
/// If fa=="", no fa comparison is performed
/// \param zip To automatically zip the output after compression for additional compression ratio
///

template<typename T>
void compress(const string & inputFileName, const string & outputFileName, const bool outOfCore, const bool timings,
              const bool displayError, const bool displayDetailedError, const bool verbose, const bool verif, const string & fa, const bool zip)
{

    std::cout << "Quantization precision: " << sizeof(T) * CHAR_BIT << " bits." << std::endl;
    Timer timer;

    if (outOfCore)
    {
        if (displayError || displayDetailedError)
            cerr << "Error can only be displayed when compressing in-core" << endl;
        if (verif)
            cout << "No verifications are done in out-of-core mode" << endl;
        if (fa != "")
            cerr << "FA computation can only be done in in-core mode" << endl;
        timer.start();
        fls::loadCompressAndSave<T>(inputFileName, outputFileName);
        if (timings)
            timer.printElapsed("Compression time");
        return;
    }

    cout << "Loading file " << inputFileName << " ... " << std::endl;
    std::vector<std::vector<Vector3f> > data;
    fls::loadTckFibers(inputFileName, data);
    cout << "Compressing fibers ..." << endl;
    timer.start();
    float ratio;
    if (verif)
        ratio = fc::computeRatioVerif(data, verbose);
    else
        ratio = fc::computeRatio(data, verbose);
    std::vector<CompressedFiber<T>> compressedFibers;
    compressedFibers.resize(data.size());
#pragma omp parallel for
    for(unsigned i = 0; i < data.size(); ++i)
        fc::compressFiber(data[i], ratio, compressedFibers[i]);
    if (timings)
        timer.printElapsed("Compression time");

    cout << "Saving file " << outputFileName << endl;
    fls::saveCompressedFibers(outputFileName, compressedFibers, ratio);


    if (zip)
    {
        cout << "Zipping file " << outputFileName << ".7z ..." << endl;
        string command = "7z a " + outputFileName + ".7z " + outputFileName;
        cout << "7z compression returned: " << system(command.c_str()) << endl;
    }

    if (fa != "")
    {
        cout << "Decompressing fibers for fa comparison ..." << endl;
        timer.stop();
        timer.start();
        std::vector<std::vector<Vector3f> > decompressedFibers;
        decompressedFibers.resize(compressedFibers.size());
#pragma omp parallel for
        for(unsigned i = 0; i < compressedFibers.size(); ++i)
            fc::decompressFiber(compressedFibers[i], ratio, decompressedFibers[i]);
        if (timings)
            timer.printElapsed("Decompression time");
        vc::compareFA(data, decompressedFibers, fa);
    }
    if (displayError)
    {
        cout << "Decompressing fibers for error computation ..." << endl;
        timer.stop();
        timer.start();
        std::vector<std::vector<Vector3f> > decompressedFibers;
        decompressedFibers.resize(compressedFibers.size());
#pragma omp parallel for
        for(unsigned i = 0; i < compressedFibers.size(); ++i)
            fc::decompressFiber(compressedFibers[i], ratio, decompressedFibers[i]);
        if (timings)
            timer.printElapsed("Decompression time");
        float maxError, meanError;
        utl::meanMaxError(data, decompressedFibers, meanError, maxError);
        cout << "Max error: " << maxError << " mm" << endl;
        cout << "Mean error: " << meanError << " mm" << endl;
    }
    if (displayDetailedError)
    {
        cout << "Decompressing fibers for detailed error computation ..." << endl;
        timer.stop();
        timer.start();
        std::vector<std::vector<Vector3f> > decompressedFibers;
        decompressedFibers.resize(compressedFibers.size());
#pragma omp parallel for
        for(unsigned i = 0; i < compressedFibers.size(); ++i)
            fc::decompressFiber(compressedFibers[i], ratio, decompressedFibers[i]);
        if (timings)
            timer.printElapsed("Decompression time");
        vector<float> minError, meanError, maxError;
        float minEndpointError, meanEndpointError, maxEndpointError;
        utl::error(data, decompressedFibers, minError, meanError, maxError, minEndpointError, meanEndpointError, maxEndpointError, 10, 40, 260);
        cout << "Max error endpoints: " << maxEndpointError << endl;
        cout << "Mean error endpoints: " << meanEndpointError << endl;
        cout << "Min error endpoints: " << minEndpointError << endl;

        fstream file;
        string outputStat =outputFileName + "-stats.txt";
        file.open(outputStat, ios_base::out);
        for (unsigned i = 0; i < maxError.size(); i++)
        {
            file << 40 + i * 10 << " ";
            file << maxError[i] << " ";
            file << meanError[i] << " ";
            file << minError[i] << endl;
        }
    }
}

///
/// \brief decompress Function to decompress a set of fibers
/// \param inputFileName Path of the input file (qfib)
/// \param outputFileName Path of the output file (tck)
/// \param outOfCore Boolean : true when outOfCore, false otherwise
/// \param timings When true, display the compression timings
///

template<typename T>
void decompress(const string & inputFileName, const string & outputFileName, const bool outOfCore, const bool timings)
{
    std::cout << "Precision of quantization used to encode this file: " << sizeof(T) * CHAR_BIT << " bits." << std::endl;
    Timer timer;

    if (outOfCore)
    {
        timer.start();
        fls::loadDecompressAndSave<T>(inputFileName, outputFileName);
        if (timings)
            timer.printElapsed("Decompression time");
        return;
    }

    cout << "Loading file " << inputFileName << endl;
    float ratio;
    std::vector<CompressedFiber<T>> compressedFibers;
    fls::loadCompressedFibers(inputFileName, compressedFibers, ratio);

    cout << "Decompressing fibers ..." << endl;
    timer.start();
    std::vector<std::vector<Vector3f> > decompressedFibers;
    decompressedFibers.resize(compressedFibers.size());
#pragma omp parallel for
    for(unsigned i = 0; i < compressedFibers.size(); ++i)
        fc::decompressFiber(compressedFibers[i], ratio, decompressedFibers[i]);
    if (timings)
        timer.printElapsed("Decompression time");

    cout << "Saving file " << outputFileName << endl;
    fls::saveTck(outputFileName, decompressedFibers);
}

///
/// \brief hasEnding Function to verify that a given file path ends with a given extension
/// \param path Path of the file to evaluate
/// \param extension Extension that the file needs to have for the function to return true
/// \return Return true if the path ends with the given extension, false otherwise
///

bool hasEnding(const string & path, const string & extension)
{
    if (path.length()<extension.length()) return false;
    return (path.compare(path.length()-extension.length(), extension.length(), extension) == 0);
}

///
/// \brief printUsage Print the simplified usage of the program
/// \param prog Name of the program
///

void printUsage(const string & prog)
{
    cerr << "Usage: " << prog << " NAME_OF_INPUT_FILE [NAME_OF_OUTPUT_FILE] [OPTIONS]\n"
                                 "See help (-h) for more informations" << endl;
}

///
/// \brief printHelp Print the help of the program
/// \param prog Name of the program
///

void printHelp(const string & prog)
{
    cout << "--------------Help-----------------" << endl;
    cout << "Usage: " <<prog << " NAME_OF_INPUT_FILE [NAME_OF_OUTPUT_FILE] [OPTIONS]" << endl;
    cout << "The input file will be either compressed or decompressed depending on its format" << endl;
    cout << "The output file will have the same name than the input file but the extension will change" << endl;
    cout << "Formats supported :" << endl;
    cout << "tck: fibers not compressed" << endl;
    cout << "qfib: fibers compressed" << endl;
    cout << "\n-------Options---------" << endl;
    cout << "-u or --out_of_core: force out-of-core, in case there is not enough RAM (default is in core)" << endl;
    cout << "-o or --octahedral: force octahedral quantification, only for compression (default is fibonacci)" << endl;
    cout << "-b or --16bits: compress using 16 bits (default is 8 bits)" << endl;
    cout << "-t or --timings: display compression and decompression times" << endl;
    cout << "-e or --error: display the error when compressing fibers" << endl;
    cout << "-d or --detailed-error: display a detailed error when compressing fibers" << endl;
    cout << "-v or --verbose: display additional details" << endl;
    cout << "-n or --no-verif: no verification of data constant stepsize for better performance" << endl;
    cout << "-f or --fa + filename: compare fa values after compression and decompression" << endl;
    cout << "-c or --compression-ratio: display compression ratio when compressing fibers" << endl;
    cout << "-h or --help: access to this help" << endl;
    cout << "Code by Corentin MERCIER and Sylvain ROUSSEAU\n" << endl;
    cout << "This source code is provided under the MIT license." << endl;
    cout << "-----------------------------------" << endl;
}

///
/// \brief main Main function, reading the arguments to perform the computation asked
/// \param argc Number of arguments
/// \param argv Char* containing the input as text
/// \return
///

int main(int argc, char* argv[])
{
//    fls::saveDifferenceAsVTK("/media/Donnees/TestTracto/SD_STREAM/tracto60kSD_STREAM1.tck","pireCasMieux.tck");
//    exit(EXIT_SUCCESS);

//    std::vector<std::vector<Vector3f>> originalFibers, decompressedFibers;
//    fls::loadTckFibers(argv[1], originalFibers);
//    fls::loadTckFibers(argv[2], decompressedFibers);
//    vc::compareFA(originalFibers, decompressedFibers, argv[3]);
//    exit(EXIT_SUCCESS);

    int c = 0;
    string inputFileName = "";
    string outputFileName = "";
    bool outOfCore = false;
    bool timings = false;
    bool displayError = false;
    bool displayDetailedError = false;
    bool bits_16 = false;
    bool verbose = false;
    bool verif = true;
    bool compressionRatio = false;
    string fa = "";

    static struct option long_options[] = {
    {"out_of_core",       no_argument,       0,  'u' },
    {"octahedral",        no_argument,       0,  'o' },
    {"16bits",            no_argument,       0,  'b' },
    {"help",              no_argument,       0,  'h' },
    {"timings",           no_argument,       0,  't' },
    {"error",             no_argument,       0,  'e' },
    {"detailed-error",    no_argument,       0,  'd' },
    {"verbose",           no_argument,       0,  'v' },
    {"no-verif",          no_argument,       0,  'n' },
    {"fa",                          1,       0,  'f' },
    {"compression-ratio", no_argument,       0,  'c' },
    {0,                             0,       0,   0  }};

    int option_index = 0;
    while ((c =  getopt_long(argc, argv, "uotedhbvnf:cj", long_options, &option_index)) != -1) {
        switch (c) {
        case 'h':
            printHelp(argv[0]);
            exit(EXIT_SUCCESS);
        case 'u':
            outOfCore = true;
            break;
        case 'o':
            quantizationMethod = QuantizationMethod::OCTAHEDRAL;
            break;
        case 'b':
            bits_16 = true;
            break;
        case 't':
            timings = true;
            break;
        case 'e':
            displayError = true;
            break;
        case 'd':
            displayDetailedError = true;
            break;
        case 'v':
            verbose = true;
            break;
        case 'n':
            verif = false;
            break;
        case 'f':
            fa = optarg;
            break;
        case 'c':
            compressionRatio = true;
            break;
        case 'j':
            cout << "Thank you Jeremie Schertzer for the 6 cents ! ;)" << endl;
            cout << "May our raspberry pi be always cool thanks to you ! 8)" << endl;
            break;
        default:
            cerr << "Character not recognized, error code: " << c << endl << endl;
            printUsage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    unsigned pos = 0;
    for (int index = optind; index < argc; index++)
    {
        switch (pos) {
        case 0:
            inputFileName = argv[index];
            break;
        case 1:
            outputFileName = argv[index];
            break;
        default:
            cerr << "Too much arguments" << endl;
            printUsage(argv[0]);
            exit(EXIT_FAILURE);
        }
        pos++;
    }

    if (inputFileName == "")
    {
        printUsage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (hasEnding(inputFileName, ".tck"))//Compression
    {
        bool zip = false;
        if (outputFileName == "")
            outputFileName = inputFileName.substr(0, inputFileName.length()-4) + ".qfib";
        if (hasEnding(outputFileName, ".7z")) {
            zip = true;
            outputFileName = outputFileName.substr(0, outputFileName.length()-3);
        }
        if (!hasEnding(outputFileName, ".qfib")){
            printUsage(argv[0]);
            exit(EXIT_FAILURE);
        }
        if (quantizationMethod == QuantizationMethod::OCTAHEDRAL)
        {
            cout << "Octahedral quantification used" << endl;
            if (bits_16)
            {
                bitcount = 16;
                compress<int16_t>(inputFileName, outputFileName, outOfCore, timings, displayError, displayDetailedError, verbose, verif, fa, zip);
            }
            else
            {
                bitcount = 8;
                compress<int8_t>(inputFileName, outputFileName, outOfCore, timings, displayError, displayDetailedError, verbose, verif, fa, zip);
            }
        }
        else
        {
            cout << "Spherical Fibonacci quantification used" << endl;
            if (bits_16)
            {
                bitcount = 16;
                compress<uint16_t>(inputFileName, outputFileName, outOfCore, timings, displayError, displayDetailedError, verbose, verif, fa, zip);
            }
            else
            {
                bitcount = 8;
                compress<uint8_t>(inputFileName, outputFileName, outOfCore, timings, displayError, displayDetailedError, verbose, verif, fa, zip);
            }
        }

        if (compressionRatio)
        {
            if (zip)
                cout << "Compression ratio: " << (1 - GetFileSize(outputFileName + ".7z") / (float)GetFileSize(inputFileName)) * 100 << "%" << endl;
            else
                cout << "Compression ratio: " << (1 - GetFileSize(outputFileName) / (float)GetFileSize(inputFileName)) * 100 << "%" << endl;
        }


    }
    else if (hasEnding(inputFileName, ".qfib") || hasEnding(inputFileName, ".qfib.7z"))//Decompression
    {
        if (hasEnding(inputFileName, ".qfib.7z")) {
            string command = "7z x " + inputFileName;
            cout << "7z decompression returned: " << system(command.c_str()) << endl;
            inputFileName = inputFileName.substr(0, inputFileName.length()-3);
        }

        if (outputFileName == "")
            outputFileName = inputFileName.substr(0, inputFileName.length()-5) + ".tck";
        if (!hasEnding(outputFileName, ".tck")){
            printUsage(argv[0]);
            exit(EXIT_FAILURE);
        }
        if (displayError || displayDetailedError)
            cerr << "Error can be displayed only whan compressing a file" << endl;
        if (fa != "")
            cerr << "No FA computation is done during decompression" << endl;
        bits_16 = fls::is16bits(inputFileName);
        if (quantizationMethod == QuantizationMethod::OCTAHEDRAL)
        {
            cout << "Octahedral quantification used" << endl;
            if (bits_16)
            {
                bitcount = 16;
                decompress<int16_t>(inputFileName, outputFileName, outOfCore, timings);
            }
            else {
                bitcount = 8;
                decompress<int8_t>(inputFileName, outputFileName, outOfCore, timings);
            }
        }
        else
        {
            cout << "Fibonacci quantification used" << endl;
            if (bits_16)
            {
                bitcount = 16;
                decompress<uint16_t>(inputFileName, outputFileName, outOfCore, timings);
            }
            else {
                bitcount = 8;
                decompress<uint8_t>(inputFileName, outputFileName, outOfCore, timings);
            }
        }
    }
    else
    {
        cerr << "Format supported :" << endl;
        cerr << "tck: fibers not compressed" << endl;
        cerr << "qfib: fibers compressed" << endl;
        cerr << "See help (-h) for more informations" << endl;
        exit(EXIT_FAILURE);
    }
    return 0;
}
