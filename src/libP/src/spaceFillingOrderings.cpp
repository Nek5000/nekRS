/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <stdlib.h>

// taken from: http://and-what-happened.blogspot.com/2011/08/fast-2d-and-3d-hilbert-curves-and.html

extern "C"
{
unsigned int MortonToHilbert2D( const unsigned int morton, const unsigned int bits )
{
  unsigned int hilbert = 0;
  unsigned int remap = 0xb4;
  unsigned int block = ( bits << 1 );
  while( block )
    {
      block -= 2;
      unsigned int mcode = ( ( morton >> block ) & 3 );
      unsigned int hcode = ( ( remap >> ( mcode << 1 ) ) & 3 );
      remap ^= ( 0x82000028 >> ( hcode << 3 ) );
      hilbert = ( ( hilbert << 2 ) + hcode );
    }
  return( hilbert );
}
unsigned int HilbertToMorton2D( const unsigned int hilbert, const unsigned int bits )
{
  unsigned int morton = 0;
  unsigned int remap = 0xb4;
  unsigned int block = ( bits << 1 );
  while( block )
    {
      block -= 2;
      unsigned int hcode = ( ( hilbert >> block ) & 3 );
      unsigned int mcode = ( ( remap >> ( hcode << 1 ) ) & 3 );
      remap ^= ( 0x330000cc >> ( hcode << 3 ) );
      morton = ( ( morton << 2 ) + mcode );
    }
  return( morton );
}
unsigned int MortonToHilbert3D( const unsigned int morton, const unsigned int bits )
{
  unsigned int hilbert = morton;
  if( bits > 1 )
    {
      unsigned int block = ( ( bits * 3 ) - 3 );
      unsigned int hcode = ( ( hilbert >> block ) & 7 );
      unsigned int mcode, shift, signs;
      shift = signs = 0;
      while( block )
        {
          block -= 3;
          hcode <<= 2;
          mcode = ( ( 0x20212021 >> hcode ) & 3 );
          shift = ( ( 0x48 >> ( 7 - shift - mcode ) ) & 3 );
          signs = ( ( signs | ( signs << 3 ) ) >> mcode );
          signs = ( ( signs ^ ( 0x53560300 >> hcode ) ) & 7 );
          mcode = ( ( hilbert >> block ) & 7 );
          hcode = mcode;
          hcode = ( ( ( hcode | ( hcode << 3 ) ) >> shift ) & 7 );
          hcode ^= signs;
          hilbert ^= ( ( mcode ^ hcode ) << block );
        }
    }
  hilbert ^= ( ( hilbert >> 1 ) & 0x92492492 );
  hilbert ^= ( ( hilbert & 0x92492492 ) >> 1 );
  return( hilbert );
}
unsigned int HilbertToMorton3D( const unsigned int hilbert, const unsigned int bits )
{
  unsigned int morton = hilbert;
  morton ^= ( ( morton & 0x92492492 ) >> 1 );
  morton ^= ( ( morton >> 1 ) & 0x92492492 );
  if( bits > 1 )
    {
      unsigned int block = ( ( bits * 3 ) - 3 );
      unsigned int hcode = ( ( morton >> block ) & 7 );
      unsigned int mcode, shift, signs;
      shift = signs = 0;
      while( block )
        {
          block -= 3;
          hcode <<= 2;
          mcode = ( ( 0x20212021 >> hcode ) & 3 );
          shift = ( ( 0x48 >> ( 4 - shift + mcode ) ) & 3 );
          signs = ( ( signs | ( signs << 3 ) ) >> mcode );
          signs = ( ( signs ^ ( 0x53560300 >> hcode ) ) & 7 );
          hcode = ( ( morton >> block ) & 7 );
          mcode = hcode;
          mcode ^= signs;
          mcode = ( ( ( mcode | ( mcode << 3 ) ) >> shift ) & 7 );
          morton ^= ( ( hcode ^ mcode ) << block );
        }
    }
  return( morton );
}
unsigned int Morton_2D_Encode_5bit( unsigned int index1, unsigned int index2 )
{ // pack 2 5-bit indices into a 10-bit Morton code
  index1 &= 0x0000001f;
  index2 &= 0x0000001f;
  index1 *= 0x01041041;
  index2 *= 0x01041041;
  index1 &= 0x10204081;
  index2 &= 0x10204081;
  index1 *= 0x00108421;
  index2 *= 0x00108421;
  index1 &= 0x15500000;
  index2 &= 0x15500000;
  return( ( index1 >> 20 ) | ( index2 >> 19 ) );
}
void Morton_2D_Decode_5bit( const unsigned int morton, unsigned int& index1, unsigned int& index2 )
{ // unpack 2 5-bit indices from a 10-bit Morton code
  unsigned int value1 = morton;
  unsigned int value2 = ( value1 >> 1 );
  value1 &= 0x00000155;
  value2 &= 0x00000155;
  value1 |= ( value1 >> 1 );
  value2 |= ( value2 >> 1 );
  value1 &= 0x00000133;
  value2 &= 0x00000133;
  value1 |= ( value1 >> 2 );
  value2 |= ( value2 >> 2 );
  value1 &= 0x0000010f;
  value2 &= 0x0000010f;
  value1 |= ( value1 >> 4 );
  value2 |= ( value2 >> 4 );
  value1 &= 0x0000001f;
  value2 &= 0x0000001f;
  index1 = value1;
  index2 = value2;
}
unsigned int Morton_2D_Encode_16bit( unsigned int index1, unsigned int index2 )
{ // pack 2 16-bit indices into a 32-bit Morton code
  index1 &= 0x0000ffff;
  index2 &= 0x0000ffff;
  index1 |= ( index1 << 8 );
  index2 |= ( index2 << 8 );
  index1 &= 0x00ff00ff;
  index2 &= 0x00ff00ff;
  index1 |= ( index1 << 4 );
  index2 |= ( index2 << 4 );
  index1 &= 0x0f0f0f0f;
  index2 &= 0x0f0f0f0f;
  index1 |= ( index1 << 2 );
  index2 |= ( index2 << 2 );
  index1 &= 0x33333333;
  index2 &= 0x33333333;
  index1 |= ( index1 << 1 );
  index2 |= ( index2 << 1 );
  index1 &= 0x55555555;
  index2 &= 0x55555555;
  return( index1 | ( index2 << 1 ) );
}
void Morton_2D_Decode_16bit( const unsigned int morton, unsigned int& index1, unsigned int& index2 )
{ // unpack 2 16-bit indices from a 32-bit Morton code
  unsigned int value1 = morton;
  unsigned int value2 = ( value1 >> 1 );
  value1 &= 0x55555555;
  value2 &= 0x55555555;
  value1 |= ( value1 >> 1 );
  value2 |= ( value2 >> 1 );
  value1 &= 0x33333333;
  value2 &= 0x33333333;
  value1 |= ( value1 >> 2 );
  value2 |= ( value2 >> 2 );
  value1 &= 0x0f0f0f0f;
  value2 &= 0x0f0f0f0f;
  value1 |= ( value1 >> 4 );
  value2 |= ( value2 >> 4 );
  value1 &= 0x00ff00ff;
  value2 &= 0x00ff00ff;
  value1 |= ( value1 >> 8 );
  value2 |= ( value2 >> 8 );
  value1 &= 0x0000ffff;
  value2 &= 0x0000ffff;
  index1 = value1;
  index2 = value2;
}
unsigned int Morton_3D_Encode_5bit( unsigned int index1, unsigned int index2, unsigned int index3 )
{ // pack 3 5-bit indices into a 15-bit Morton code
  index1 &= 0x0000001f;
  index2 &= 0x0000001f;
  index3 &= 0x0000001f;
  index1 *= 0x01041041;
  index2 *= 0x01041041;
  index3 *= 0x01041041;
  index1 &= 0x10204081;
  index2 &= 0x10204081;
  index3 &= 0x10204081;
  index1 *= 0x00011111;
  index2 *= 0x00011111;
  index3 *= 0x00011111;
  index1 &= 0x12490000;
  index2 &= 0x12490000;
  index3 &= 0x12490000;
  return( ( index1 >> 16 ) | ( index2 >> 15 ) | ( index3 >> 14 ) );
}
void Morton_3D_Decode_5bit( const unsigned int morton,
			    unsigned int& index1, unsigned int& index2, unsigned int& index3 )
{ // unpack 3 5-bit indices from a 15-bit Morton code
  unsigned int value1 = morton;
  unsigned int value2 = ( value1 >> 1 );
  unsigned int value3 = ( value1 >> 2 );
  value1 &= 0x00001249;
  value2 &= 0x00001249;
  value3 &= 0x00001249;
  value1 |= ( value1 >> 2 );
  value2 |= ( value2 >> 2 );
  value3 |= ( value3 >> 2 );
  value1 &= 0x000010c3;
  value2 &= 0x000010c3;
  value3 &= 0x000010c3;
  value1 |= ( value1 >> 4 );
  value2 |= ( value2 >> 4 );
  value3 |= ( value3 >> 4 );
  value1 &= 0x0000100f;
  value2 &= 0x0000100f;
  value3 &= 0x0000100f;
  value1 |= ( value1 >> 8 );
  value2 |= ( value2 >> 8 );
  value3 |= ( value3 >> 8 );
  value1 &= 0x0000001f;
  value2 &= 0x0000001f;
  value3 &= 0x0000001f;
  index1 = value1;
  index2 = value2;
  index3 = value3;
}
unsigned int Morton_3D_Encode_10bit( unsigned int index1, unsigned int index2, unsigned int index3 )
{ // pack 3 10-bit indices into a 30-bit Morton code
  index1 &= 0x000003ff;
  index2 &= 0x000003ff;
  index3 &= 0x000003ff;
  index1 |= ( index1 << 16 );
  index2 |= ( index2 << 16 );
  index3 |= ( index3 << 16 );
  index1 &= 0x030000ff;
  index2 &= 0x030000ff;
  index3 &= 0x030000ff;
  index1 |= ( index1 << 8 );
  index2 |= ( index2 << 8 );
  index3 |= ( index3 << 8 );
  index1 &= 0x0300f00f;
  index2 &= 0x0300f00f;
  index3 &= 0x0300f00f;
  index1 |= ( index1 << 4 );
  index2 |= ( index2 << 4 );
  index3 |= ( index3 << 4 );
  index1 &= 0x030c30c3;
  index2 &= 0x030c30c3;
  index3 &= 0x030c30c3;
  index1 |= ( index1 << 2 );
  index2 |= ( index2 << 2 );
  index3 |= ( index3 << 2 );
  index1 &= 0x09249249;
  index2 &= 0x09249249;
  index3 &= 0x09249249;
  return( index1 | ( index2 << 1 ) | ( index3 << 2 ) );
}
void Morton_3D_Decode_10bit( const unsigned int morton,
			     unsigned int& index1, unsigned int& index2, unsigned int& index3 )
{ // unpack 3 10-bit indices from a 30-bit Morton code
  unsigned int value1 = morton;
  unsigned int value2 = ( value1 >> 1 );
  unsigned int value3 = ( value1 >> 2 );
  value1 &= 0x09249249;
  value2 &= 0x09249249;
  value3 &= 0x09249249;
  value1 |= ( value1 >> 2 );
  value2 |= ( value2 >> 2 );
  value3 |= ( value3 >> 2 );
  value1 &= 0x030c30c3;
  value2 &= 0x030c30c3;
  value3 &= 0x030c30c3;
  value1 |= ( value1 >> 4 );
  value2 |= ( value2 >> 4 );
  value3 |= ( value3 >> 4 );
  value1 &= 0x0300f00f;
  value2 &= 0x0300f00f;
  value3 &= 0x0300f00f;
  value1 |= ( value1 >> 8 );
  value2 |= ( value2 >> 8 );
  value3 |= ( value3 >> 8 );
  value1 &= 0x030000ff;
  value2 &= 0x030000ff;
  value3 &= 0x030000ff;
  value1 |= ( value1 >> 16 );
  value2 |= ( value2 >> 16 );
  value3 |= ( value3 >> 16 );
  value1 &= 0x000003ff;
  value2 &= 0x000003ff;
  value3 &= 0x000003ff;
  index1 = value1;
  index2 = value2;
  index3 = value3;
}



unsigned int hilbert2D(unsigned int index1, unsigned int index2){

  unsigned int morton = Morton_2D_Encode_16bit(index1,index2);

  return MortonToHilbert2D(morton, 16);
}

unsigned int morton2D(unsigned int index1, unsigned int index2){

  return Morton_2D_Encode_16bit(index1,index2);
}



}
