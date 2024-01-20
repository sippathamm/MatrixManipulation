#pragma once

#ifndef NO_EXTRAS
#include <ostream>
#endif

#include <Vector.h>

#include <cassert>

/*
 * Matrix transformation performs in row-major.
 * Converting between row-major and column-major is taking transpose of matrix.
 *
 * Row-major:
 * |-------------|---|  |-----------|---|
 * |             | 0 |  |           | 0 |
 * |   Rotation  | 0 |  |   (3x3)   | 0 |
 * |             | 0 |  |           | 0 |
 * |-------------|---|  |-----------|---|
 * | Translation | 1 |  |   (1x3)   | 1 |
 * |-------------|---|  |-----------|---|
 *
 * Column-major:
 * |-------------|-------------|    |-----------|-------|
 * |   Rotation  | Translation |    |   (3x3)   | (3x1) |
 * |-------------|-------------|    |-----------|-------|
 * |   0  0  0   |      1      |    |  0  0  0  |   1   |
 * |-------------|-------------|    |-----------|-------|
 */

namespace SPTH
{
    typedef struct A2DMatrix
    {
        inline A2DMatrix ()
        {
            M11 = M22 = 1.0f;
            M12 = M21 = 0.0f;
        }

        inline A2DMatrix (float GivenM11, float GivenM12,
                          float GivenM21, float GivenM22)
        {
            M11 = GivenM11; M12 = GivenM12;
            M21 = GivenM21; M22 = GivenM22;
        }

        union
        {
            struct
            {
                float M11, M12,
                      M21, M22;
            };
            float MatrixArray[4] {};
        };

        inline float* operator [] (int MatrixIndex);
        A2DMatrix operator * (const A2DMatrix &AnotherMatrix) const;
        A2DMatrix operator * (float Constant) const;
        A2DMatrix Transpose () const;
        A2DMatrix Minor () const;
        A2DMatrix Cofactor () const;
        float Determinant () const;
        A2DMatrix Adjugate () const;
        A2DMatrix Inverse () const;
    } A2DMatrix;

    typedef struct A3DMatrix
    {
        inline A3DMatrix ()
        {
            M11 = M22 = M33 = 1.0f;
            M12 = M13 = M21 = 0.0f;
            M23 = M31 = M32 = 0.0f;
        }

        inline A3DMatrix (float GivenM11, float GivenM12, float GivenM13,
                          float GivenM21, float GivenM22, float GivenM23,
                          float GivenM31, float GivenM32, float GivenM33)
        {
            M11 = GivenM11; M12 = GivenM12; M13 = GivenM13;
            M21 = GivenM21; M22 = GivenM22; M23 = GivenM23;
            M31 = GivenM31; M32 = GivenM32; M33 = GivenM33;
        }

        union
        {
            struct
            {
                float M11, M12, M13,
                      M21, M22, M23,
                      M31, M32, M33;
            };
            float MatrixArray[9] {};
        };

        inline float* operator [] (int MatrixIndex);
        A3DMatrix operator * (const A3DMatrix &AnotherMatrix) const;
        A3DMatrix operator * (float Constant) const;
        A3DMatrix Transpose () const;
        A2DMatrix Cut (int Row, int Column) const;
        A3DMatrix Minor () const;
        A3DMatrix Cofactor () const;
        float Determinant () const;
        A3DMatrix Adjugate () const;
        A3DMatrix Inverse () const;
        A3DVector MultiplyVector (const A3DVector &Vector) const;
    } A3DMatrix;

    typedef struct A4DMatrix
    {
        inline A4DMatrix ()
        {
            M11 = M22 = M33 = M44 = 1.0f;
            M12 = M13 = M14 = M21 = 0.0f;
            M23 = M31 = M24 = M32 = 0.0f;
            M34 = M41 = M42 = M43 = 0.0f;
        }

        inline A4DMatrix (float GivenM11, float GivenM12, float GivenM13, float GivenM14,
                          float GivenM21, float GivenM22, float GivenM23, float GivenM24,
                          float GivenM31, float GivenM32, float GivenM33, float GivenM34,
                          float GivenM41, float GivenM42, float GivenM43, float GivenM44)
        {
            M11 = GivenM11; M12 = GivenM12; M13 = GivenM13; M14 = GivenM14;
            M21 = GivenM21; M22 = GivenM22; M23 = GivenM23; M24 = GivenM24;
            M31 = GivenM31; M32 = GivenM32; M33 = GivenM33; M34 = GivenM34;
            M41 = GivenM41; M42 = GivenM42; M43 = GivenM43; M44 = GivenM44;
        }

        union
        {
            struct
            {
                float M11, M12, M13, M14,
                      M21, M22, M23, M24,
                      M31, M32, M33, M34,
                      M41, M42, M43, M44;
            };
            float MatrixArray[16] {};
        };

        inline float* operator [] (int MatrixIndex)
        {
            assert (MatrixIndex >= 0 && MatrixIndex < 4);
            return &(MatrixArray[MatrixIndex * 4]);
        }
        A4DMatrix operator * (const A4DMatrix &AnotherMatrix) const;
        A4DMatrix operator * (float Constant) const;
        A4DMatrix Transpose () const;
        A3DMatrix Cut (int Row, int Column) const;
        A4DMatrix Minor () const;
        A4DMatrix Cofactor () const;
        float Determinant () const;
        A4DMatrix Adjugate () const;
        A4DMatrix Inverse () const;
        A3DVector GetTranslation () const;
        A3DVector GetScale () const;
        A3DVector MultiplyPoint (const A3DVector &Point) const;
        A3DVector MultiplyVector (const A3DVector &Vector) const;
    } A4DMatrix;

    void Transpose (float *ResultMatrixArray, const float *Matrix, int Row, int Column);
    void Multiply (float *ResultMatrixArray, const float *MatrixA, int MatrixARow, int MatrixAColumn, const float *MatrixB, int MatrixBRow, int MatrixBColumn);
    void Cofactor (float *ResultMatrixArray, const float *Matrix, int Row, int Column);

    A3DMatrix Rotation3D (float XAngle, float YAngle, float ZAngle);
    A3DMatrix ZRotation3D (float ZAngle);
    A3DMatrix XRotation3D (float XAngle);
    A3DMatrix YRotation3D (float YAngle);
    A3DMatrix AxisAngle3D (const A3DVector &Axis, float Angle);

    A4DMatrix Translation (float X, float Y, float Z);
    A4DMatrix Translation (const A3DVector &TranslationVector);
    A4DMatrix FromMatrix3D (const A3DMatrix &Matrix);
    A4DMatrix Scale (float X, float Y, float Z);
    A4DMatrix Scale (const A3DVector &ScaleVector);
    A4DMatrix Rotation4D (float XAngle, float YAngle, float ZAngle);
    A4DMatrix ZRotation4D (float ZAngle);
    A4DMatrix XRotation4D (float XAngle);
    A4DMatrix YRotation4D (float YAngle);
    A4DMatrix AxisAngle4D (const A3DVector &Axis, float Angle);
    A4DMatrix Transformation (const A3DVector &ScaleVector, const A3DVector &EulerRotationVector, const A3DVector &TranslationVector);
    A4DMatrix Transformation (const A3DVector &ScaleVector, const A3DVector &AxisRotationVector, float Angle, const A3DVector &TranslationVector);
    A4DMatrix LookAt (const A3DVector &CameraPosition, const A3DVector &TargetPosition, const A3DVector &UpwardVector);
    A4DMatrix Perspective (float FOV, float Aspect, float ZNear, float ZFar);
    A4DMatrix Orthographic (float Left, float Right, float Bottom, float Top, float ZNear, float ZFar);
}

#ifndef NO_EXTRAS
std::ostream& operator << (std::ostream &OS, const SPTH::A2DMatrix &MatrixA);
std::ostream& operator << (std::ostream &OS, const SPTH::A3DMatrix &MatrixA);
std::ostream& operator << (std::ostream &OS, const SPTH::A4DMatrix &MatrixA);
#endif