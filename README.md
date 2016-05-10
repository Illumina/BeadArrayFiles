# IlluminaBeadArrayFiles
Library to parse file formats related to Illumina bead arrays

## Build instructions
Need to provide instructions for building package with easy install

## GTC File Format
The file format is composed of a “Table of Contents” at the head of the file. The table of contents has an arbitrary number of tables of content
(TOC) entries. This allows for forwards and backwards compatibility as new fields are added to the file format or old fields are removed. A
table-of-content entry is 6 bytes long and is made up of two parts: an ID and a file OFFSET. The ID is a two-byte short integer that identifies what
the variable is, and the file OFFSET is a four-byte integer that identifies the file offset (relative to the start of the file) where the data for the
variable resides. Table 1 describes the format for the file header. Table 2 describes the format for a TOC entry. Table 3 describes the IDs
recognized in the gtc file specification. **Note that all multi-byte variables are stored with lowest order byte first.**

####Table 1: File format header
<table>
   <tbody>
      <tr>
         <th>
            <p>Parameter </p>
         </th>
         <th>
            <p>File Offset</p>
         </th>
         <th>
            <p>Type</p>
         </th>
         <th>
            <p>Description</p>
         </th>
      </tr>
      <tr>
         <td>
            <p>Identifier</p>
         </td>
         <td>
            <p>0</p>
         </td>
         <td>
            <p>char[3]</p>
         </td>
         <td colspan="2">
            <p>&lsquo;gtc&rsquo;</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>File Version</p>
         </td>
         <td>
            <p>3</p>
         </td>
         <td>
            <p>byte</p>
         </td>
         <td colspan="2">
            <p>File version (5 is the version corresponding to this document)</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>Num Entries</p>
         </td>
         <td>
            <p>4</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td colspan="2">
            <p>Number of table-of-content entries (M)</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>Entry 0</p>
         </td>
         <td>
            <p>8</p>
         </td>
         <td>
            <p>TOC Entry</p>
         </td>
         <td colspan="2">
            <p>Table-of-content entry 0</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>&hellip;</p>
         </td>
         <td>
            <p> </p>
         </td>
         <td></td>
         <td colspan="2">
            <p>Table-of-content entries</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>Entry M-1</p>
         </td>
         <td>
            <p>8 + (M-1) * 6</p>
         </td>
         <td>
            <p>TOC Entry</p>
         </td>
         <td colspan="2">
            <p>Table-of-content entry M-1</p>
         </td>
      </tr>
   </tbody>
</table>
<br>

####Table 2: TOC entry
<table>
   <thead>
      <tr>
         <th tabindex="0" data-column="0">
            <div>
               <p>Parameter </p>
            </div>
         </th>
         <th tabindex="0" data-column="1">
            <div>
               <p>Type</p>
            </div>
         </th>
         <th tabindex="0" data-column="2">
            <div>
               <p>Description</p>
            </div>
         </th>
      </tr>
   </thead>
   <tbody>
      <tr>
         <td>
            <p>ID</p>
         </td>
         <td>
            <p>short</p>
         </td>
         <td>
            <p>TOC variable ID</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>OFFSET</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td>
            <p>Offset into the file (relative to start of file)</p>
         </td>
      </tr>
   </tbody>
</table>
<br>

####Table 3: Description of TOC variable IDs
<table>
   <thead>
      <tr>
         <th tabindex="0" data-column="0">
            <div>
               <p>ID</p>
            </div>
         </th>
         <th tabindex="0" data-column="1">
            <div>
               <p>Name</p>
            </div>
         </th>
         <th tabindex="0" data-column="2">
            <div>
               <p>Type</p>
            </div>
         </th>
         <th tabindex="0" data-column="3">
            <div>
               <p>Description</p>
            </div>
         </th>
      </tr>
   </thead>
   <tbody>
      <tr>
         <td>
            <p>1</p>
         </td>
         <td>
            <p>NumSNPs</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td>
            <p>Number of SNPs</p>
         </td>
      </tr>
      <tr>
         <td colspan="1">2</td>
         <td colspan="1">Ploidy</td>
         <td colspan="1">int32</td>
         <td colspan="1">Ploidy of species</td>
      </tr>
      <tr>
         <td colspan="1">3</td>
         <td colspan="1">Ploidy Type</td>
         <td colspan="1">int32</td>
         <td colspan="1">
            <div>
               <table>
                  <thead>
                     <tr>
                        <th tabindex="0" data-column="0">
                           <div>Value</div>
                        </th>
                        <th tabindex="0" data-column="1">
                           <div>Description</div>
                        </th>
                     </tr>
                  </thead>
                  <tbody>
                     <tr>
                        <td>1</td>
                        <td>Diploid</td>
                     </tr>
                     <tr>
                        <td>2</td>
                        <td>Autopolyploid</td>
                     </tr>
                     <tr>
                        <td colspan="1">3</td>
                        <td colspan="1">Allopolyploid</td>
                     </tr>
                  </tbody>
               </table>
            </div>
         </td>
      </tr>
      <tr>
         <td>
            <p>10</p>
         </td>
         <td>
            <p>Sample Name</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The sample name</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>11</p>
         </td>
         <td>
            <p>Sample Plate</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The sample plate</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>12</p>
         </td>
         <td>
            <p>Sample Well</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The sample well</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>100</p>
         </td>
         <td>
            <p>Cluster File</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The name of the cluster file</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>101</p>
         </td>
         <td>
            <p>SNP Manifest</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The name of the SNP manifest</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>200</p>
         </td>
         <td>
            <p>Imaging Date</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The imaging date</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>201</p>
         </td>
         <td>
            <p>AutoCall Date</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The AutoCall processing date</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>300</p>
         </td>
         <td>
            <p>AutoCall Version</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The AutoCall version</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>400</p>
         </td>
         <td>
            <p>Normalization Transformations</p>
         </td>
         <td>
            <p>XForm[]</p>
         </td>
         <td>
            <p>The array of normalization transformation</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>500</p>
         </td>
         <td>
            <p>Raw Control X</p>
         </td>
         <td>
            <p>ushort[]</p>
         </td>
         <td>
            <p>Raw control green intensities. The length of the vector will be equal to the (sections/sample) * (number of controls). Each control intensity will then be repeated for each section per sample.</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>501</p>
         </td>
         <td>
            <p>Raw Control Y</p>
         </td>
         <td>
            <p>ushort[]</p>
         </td>
         <td>
            <p>Raw control red intensities. The length of the vector will be equal to the (sections/sample) * (number of controls). Each control intensity will then be repeated for each section per sample.</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1000</p>
         </td>
         <td>
            <p>Raw X intensity values</p>
         </td>
         <td>
            <p>ushort[]</p>
         </td>
         <td>
            <p>The array of raw green intensities for every SNP</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1001</p>
         </td>
         <td>
            <p>Raw Y intensity values</p>
         </td>
         <td>
            <p>ushort[]</p>
         </td>
         <td>
            <p>The array of raw red intensities for every SNP</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1002</p>
         </td>
         <td>
            <p>Genotypes</p>
         </td>
         <td>
            <p>byte[]</p>
         </td>
         <td>
            <p>The array of genotypes (see table 4: genotype mapping below) for every SNP</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1003</p>
         </td>
         <td>
            <p>BaseCalls</p>
         </td>
         <td>
            <p>BaseCall[]</p>
         </td>
         <td>
            <p>The array of BaseCalls with respect to the TOP strand for every SNP</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1004</p>
         </td>
         <td>
            <p>Genotype Scores</p>
         </td>
         <td>
            <p>float[]</p>
         </td>
         <td>
            <p>The array of GenCall scores for every SNP</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1005</p>
         </td>
         <td>
            <p><em>Scanner Data</em></p>
         </td>
         <td>
            <p></p>
         </td>
         <td>
            <p><em>Information about the scanner</em></p>
         </td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>Scanner Name</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>The name of the scanner</p>
         </td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>Pmt Green</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td></td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>Pmt Red</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td></td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>Scanner Version</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>Version of the scanner software used</p>
         </td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>Imaging User</p>
         </td>
         <td>
            <p>string</p>
         </td>
         <td>
            <p>Name of the scanner user</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1006</p>
         </td>
         <td>
            <p>Call Rate</p>
         </td>
         <td>
            <p>float</p>
         </td>
         <td>
            <p>Calculated call rate of the sample</p>
         </td>
      </tr>
      <tr>
         <td>
            <p>1007</p>
         </td>
         <td>
            <p>Estimated Gender</p>
         </td>
         <td>
            <p>char</p>
         </td>
         <td>
            <p>&lsquo;M&rsquo;-Male, &lsquo;F&rsquo;-Female, &lsquo;U&rsquo;-Unknown</p>
         </td>
      </tr>
      <tr>
         <td colspan="1">1008</td>
         <td colspan="1">LogR Dev</td>
         <td colspan="1">float</td>
         <td></td>
      </tr>
      <tr>
         <td>
            <p>1009</p>
         </td>
         <td>
            <p>p10GC</p>
         </td>
         <td>
            <p>float</p>
         </td>
         <td>
            <p>The 10th percentile of genotype call scores for this sample. In the calculation of this metrics, scores from intensity only loci or scores beneath the genotype call threshold are <strong>not</strong> considered.</p>
         </td>
      </tr>
      <tr>
         <td colspan="1">1010</td>
         <td colspan="1">DX</td>
         <td colspan="1">int32</td>
         <td colspan="1">DX flag</td>
      </tr>
      <tr>
         <td>
            <p>1011</p>
         </td>
         <td>
            <p><em>Extended sample data</em></p>
         </td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">P50GC</td>
         <td colspan="1">float</td>
         <td colspan="1">The 50th percentile of genotype call scores for this sample. In the calculation of this metrics, scores from intensity only loci or scores beneath the genotype call threshold are <strong>not</strong> considered.</td>
         <td></td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>NumCalls</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td>
            <p>Number of valid calls</p>
         </td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>NumNoCalls</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td>
            <p>Number of invalid calls</p>
         </td>
      </tr>
      <tr>
         <td></td>
         <td>
            <p>Num Intensity Only</p>
         </td>
         <td>
            <p>int32</p>
         </td>
         <td>
            <p>Number of loci that are &ldquo;Intensity Only&rdquo; or "Zeroed"</p>
         </td>
      </tr>
      <tr>
         <td colspan="1">1012</td>
         <td colspan="1">B Allele Frequencies</td>
         <td colspan="1">float[]</td>
         <td colspan="1">B allele frequencies across loci</td>
      </tr>
      <tr>
         <td colspan="1">1013</td>
         <td colspan="1">LogR Ratios</td>
         <td colspan="1">float[]</td>
         <td colspan="1">LogR ratios across loci</td>
      </tr>
      <tr>
         <td colspan="1">1014</td>
         <td colspan="1"><em>Intensity percentiles (X)</em></td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">P05 X</td>
         <td colspan="1">ushort</td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">P50 X</td>
         <td colspan="1">ushort</td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">P95 X</td>
         <td colspan="1">ushort</td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">1015</td>
         <td colspan="1"><em>Intensity percentiles (Y)</em></td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">P05 Y</td>
         <td colspan="1">ushort</td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">P50 Y</td>
         <td colspan="1">ushort</td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">P95 Y</td>
         <td colspan="1">ushort</td>
         <td></td>
         <td></td>
      </tr>
      <tr>
         <td colspan="1">1016</td>
         <td colspan="1">Sentrix ID</td>
         <td colspan="1">string</td>
         <td colspan="1">The Sentrix barcode.</td>
      </tr>
   </tbody>
</table>
<br>
####Table 4: Genotype mapping table
<div>
   <table>
      <thead>
         <tr>
            <th tabindex="0" data-column="0">
               <div>Genotype</div>
            </th>
            <th tabindex="0" data-column="1">
               <div>Byte value</div>
            </th>
         </tr>
      </thead>
      <tbody>
         <tr>
            <td>NC</td>
            <td>0</td>
         </tr>
         <tr>
            <td>AA</td>
            <td>1</td>
         </tr>
         <tr>
            <td>AB</td>
            <td>2</td>
         </tr>
         <tr>
            <td>BB</td>
            <td>3</td>
         </tr>
         <tr>
            <td>A</td>
            <td>4</td>
         </tr>
         <tr>
            <td>B</td>
            <td>5</td>
         </tr>
         <tr>
            <td>AAA</td>
            <td>6</td>
         </tr>
         <tr>
            <td colspan="1">AAB</td>
            <td colspan="1">7</td>
         </tr>
         <tr>
            <td colspan="1">ABB</td>
            <td colspan="1">8</td>
         </tr>
         <tr>
            <td colspan="1">BBB</td>
            <td colspan="1">9</td>
         </tr>
         <tr>
            <td colspan="1">AAAA</td>
            <td colspan="1">10</td>
         </tr>
         <tr>
            <td colspan="1">AAAB</td>
            <td colspan="1">11</td>
         </tr>
         <tr>
            <td colspan="1">AABB</td>
            <td colspan="1">12</td>
         </tr>
         <tr>
            <td colspan="1">ABBB</td>
            <td colspan="1">13</td>
         </tr>
         <tr>
            <td colspan="1">BBBB</td>
            <td colspan="1">14</td>
         </tr>
         <tr>
            <td colspan="1">...</td>
            <td colspan="1">...</td>
         </tr>
         <tr>
            <td colspan="1">BBBBBBBB</td>
            <td colspan="1">44</td>
         </tr>
      </tbody>
   </table>
</div>

If the ID in the TOC entry corresponds to an int variable type (NumSNPs or Ploidy, in our case), then the OFFSET in the TOC entry is the actual
value of the unsigned integer and not a file offset. If the ID in the TOC entry corresponds to a string variable type, then the first byte at the file
location specified by OFFSET is the length of the string (L) and the L bytes of the string follow as single byte characters. If the ID in the TOC
corresponds to an array, then the first four bytes at the location specified by OFFSET are an integer corresponding to the length of the array (N).
Each element in the array follows. Arrays of type ushort have N ushorts. Arrays of type byte have N bytes. Arrays of type float have N floats.
The other two types require a bit of explanation. 

###BaseCall
A base call is a multiple-character word. The number of characters in the basecall is always two, regardless of the ploidy. The characters are ‘A’,
‘C’, ‘G’, ‘T’, or ‘-‘ for a no-call. For a diploid genotype, the BaseCall will be the nucleotide genotype (relative to the top strand). For a polyploid
genotype, the BaseCall will still be exactly two characters. In this case, the software will report the nucleotide genotype for the A and B alleles on
the top strand (in that order). In combination with the AB genotypes, the client will be able to reconstruct the full nucleotide genotype on the top
strand. 

###XForm
An XForm is a 52 byte structure made up of 1 4-byte integer followed by 12 4-byte floating point numbers, summarized in Table 4.

####Table 5: Normalization Transformation structure
<div>
   <table>
      <thead>
         <tr>
            <th tabindex="0" colspan="1" data-column="0">
               <div>Offset</div>
            </th>
            <th tabindex="0" colspan="1" data-column="1">
               <div>Column Name</div>
            </th>
            <th tabindex="0" colspan="1" data-column="2">
               <div>Type</div>
            </th>
         </tr>
      </thead>
      <tbody>
         <tr>
            <td>0</td>
            <td>version</td>
            <td>int</td>
         </tr>
         <tr>
            <td>4</td>
            <td>offset_x</td>
            <td>float</td>
         </tr>
         <tr>
            <td>8</td>
            <td>offset_y</td>
            <td>float</td>
         </tr>
         <tr>
            <td>12</td>
            <td>scale_x</td>
            <td>float</td>
         </tr>
         <tr>
            <td>16</td>
            <td>scale_y</td>
            <td>float</td>
         </tr>
         <tr>
            <td>20</td>
            <td>shear</td>
            <td>float</td>
         </tr>
         <tr>
            <td>24</td>
            <td>theta</td>
            <td>float</td>
         </tr>
         <tr>
            <td>28</td>
            <td>reserved</td>
            <td>float</td>
         </tr>
         <tr>
            <td>32</td>
            <td>reserved</td>
            <td>float</td>
         </tr>
         <tr>
            <td>36</td>
            <td>reserved</td>
            <td>float</td>
         </tr>
         <tr>
            <td>40</td>
            <td>reserved</td>
            <td>float</td>
         </tr>
         <tr>
            <td>44</td>
            <td>reserved</td>
            <td>float</td>
         </tr>
         <tr>
            <td>48</td>
            <td>reserved</td>
            <td>float</td>
         </tr>
      </tbody>
   </table>
</div>

To go from raw coordinates (xraw, yraw) to normalized coordinates (xn, yn), perform the following operations:

<div>
   <p>
      tempx = xraw - offset_x
      <br />
      tempy = yraw - offset_y
   </p>
   <p>
      tempx2 = cos(theta) * tempx + sin(theta) * tempy
      <br />
      tempy2 = -sin(theta) * tempx + cos(theta) * tempy
   </p>
   <p>
      tempx3 = tempx2 - shear * tempy2
      <br />
      tempy3 = tempy2
   </p>
   <p>
      xn = tempx3 / scale_x
      <br />
      yn = tempy3 / scale_y
   </p>
</div>

The reason that there is an array of normalization transformations is that a given normalization transformation applies to a subset of the SNPs in a
product. You can use the NormID parameter from the map file to identify which normalization transformation to apply to every SNP. 
