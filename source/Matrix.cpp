#include <Matrix.h>

#include <cmath>
#include <cfloat>

#define COMPARISON(A, B) (fabsf((A) - (B)) <= FLT_EPSILON * fmaxf(1.0f, fmaxf(fabsf(A), fabsf(B))))

inline float* SPTH::A2DMatrix::operator [] (int MatrixIndex)
{
    assert (MatrixIndex >= 0 && MatrixIndex < 2);
    
    return &(MatrixArray[MatrixIndex * 2]);
}

SPTH::A2DMatrix SPTH::A2DMatrix::operator * (const SPTH::A2DMatrix &AnotherMatrix) const
{
    A2DMatrix ResultMatrix;
    SPTH::Multiply(ResultMatrix.MatrixArray, this->MatrixArray, 2, 2, AnotherMatrix.MatrixArray, 2, 2);
    
    return ResultMatrix;
}

SPTH::A2DMatrix SPTH::A2DMatrix::operator * (float Constant) const
{
    A2DMatrix ResultMatrix;
    
    for (int Index = 0; Index < 4; ++Index)
    {
        ResultMatrix.MatrixArray[Index] = this->MatrixArray[Index] * Constant;
    }
    
    return ResultMatrix;
}

SPTH::A2DMatrix SPTH::A2DMatrix::Transpose () const
{
    A2DMatrix ResultMatrix;
    SPTH::Transpose(ResultMatrix.MatrixArray, this->MatrixArray, 2, 2);
    
    return ResultMatrix;
}

SPTH::A2DMatrix SPTH::A2DMatrix::Minor () const
{
    return {this->M22, this->M21, this->M12, this->M11};
}

SPTH::A2DMatrix SPTH::A2DMatrix::Cofactor () const
{
    A2DMatrix ResultMatrix;
    SPTH::Cofactor(ResultMatrix.MatrixArray, this->Minor().MatrixArray, 2, 2);
    
    return ResultMatrix;
}

float SPTH::A2DMatrix::Determinant () const
{
    return this->M11 * this->M22 - this->M21 * this->M12;
}

SPTH::A2DMatrix SPTH::A2DMatrix::Adjugate () const
{
    return this->Cofactor().Transpose();
}

SPTH::A2DMatrix SPTH::A2DMatrix::Inverse () const
{
    float Determinant = this->Determinant();
    
    if (COMPARISON(Determinant, 0.0f))
    {
        return *this;
    }
    
    return this->Adjugate() * (1.0f / Determinant);
}

inline float* SPTH::A3DMatrix::operator [] (int MatrixIndex)
{
    assert (MatrixIndex >= 0 && MatrixIndex < 3);
    
    return &(MatrixArray[MatrixIndex * 3]);
}

SPTH::A3DMatrix SPTH::A3DMatrix::operator * (const SPTH::A3DMatrix &AnotherMatrix) const
{
    A3DMatrix ResultMatrix;
    SPTH::Multiply(ResultMatrix.MatrixArray, this->MatrixArray, 3, 3, AnotherMatrix.MatrixArray, 3, 3);
    
    return ResultMatrix;
}

SPTH::A3DMatrix SPTH::A3DMatrix::operator * (float Constant) const
{
    A3DMatrix ResultMatrix;
    
    for (int Index = 0; Index < 9; ++Index)
    {
        ResultMatrix.MatrixArray[Index] = this->MatrixArray[Index] * Constant;
    }
    
    return ResultMatrix;
}

SPTH::A3DMatrix SPTH::A3DMatrix::Transpose () const
{
    A3DMatrix ResultMatrix;
    SPTH::Transpose(ResultMatrix.MatrixArray, this->MatrixArray, 3, 3);
    
    return ResultMatrix;
}

SPTH::A2DMatrix SPTH::A3DMatrix::Cut (int Row, int Column) const
{
    A2DMatrix ResultMatrix;
    int Index = 0;
    
    for (int RowIndex = 0; RowIndex < 3; ++RowIndex)
    {
        for (int ColumnIndex = 0; ColumnIndex < 3; ++ColumnIndex)
        {
            if (RowIndex == Row || ColumnIndex == Column)
            {
                continue;
            }
            ResultMatrix.MatrixArray[Index++] = this->MatrixArray[3 * RowIndex + ColumnIndex];
        }
    }
    
    return ResultMatrix;
}

SPTH::A3DMatrix SPTH::A3DMatrix::Minor () const
{
    A3DMatrix ResultMatrix;
    
    for (int RowIndex = 0; RowIndex < 3; ++RowIndex)
    {
        for (int ColumnIndex = 0; ColumnIndex < 3; ++ColumnIndex)
        {
            ResultMatrix[RowIndex][ColumnIndex] = this->Cut(RowIndex, ColumnIndex).Determinant();
        }
    }
    
    return ResultMatrix;
}

SPTH::A3DMatrix SPTH::A3DMatrix::Cofactor () const
{
    A3DMatrix ResultMatrix;
    SPTH::Cofactor(ResultMatrix.MatrixArray, this->Minor().MatrixArray, 3, 3);
    
    return ResultMatrix;
}

float SPTH::A3DMatrix::Determinant () const
{
    float ResultDeterminant = 0.0f;
    A3DMatrix ResultCofactorMatrix = this->Cofactor();
    
    for (int Index = 0; Index < 3; ++Index)
    {
        ResultDeterminant += this->MatrixArray[Index] * ResultCofactorMatrix[0][Index];
    }
    
    return ResultDeterminant;
}

SPTH::A3DMatrix SPTH::A3DMatrix::Adjugate () const
{
    return this->Cofactor().Transpose();
}

SPTH::A3DMatrix SPTH::A3DMatrix::Inverse () const
{
    float ResultDeterminant = this->Determinant();
    
    if (COMPARISON(ResultDeterminant, 0.0f))
    {
        return *this;
    }
    
    return this->Adjugate() * (1.0f / ResultDeterminant);
}

SPTH::A3DVector SPTH::A3DMatrix::MultiplyVector (const A3DVector &Vector) const
{
    A3DVector ResultVector;
    ResultVector.X = Vector.Dot({this->M11, this->M21, this->M31});
    ResultVector.Y = Vector.Dot({this->M12, this->M22, this->M32});
    ResultVector.Z = Vector.Dot({this->M13, this->M23, this->M33});
    
    return ResultVector;
}

SPTH::A4DMatrix SPTH::A4DMatrix::operator * (const SPTH::A4DMatrix &AnotherMatrix) const
{
    A4DMatrix ResultMatrix;
    SPTH::Multiply(ResultMatrix.MatrixArray, this->MatrixArray, 4, 4, AnotherMatrix.MatrixArray, 4, 4);
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::A4DMatrix::operator * (float Constant) const
{
    A4DMatrix ResultMatrix;
    
    for (int Index = 0; Index < 16; ++Index)
    {
        ResultMatrix.MatrixArray[Index] = this->MatrixArray[Index] * Constant;
    }
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::A4DMatrix::Transpose () const
{
    A4DMatrix ResultMatrix;
    SPTH::Transpose(ResultMatrix.MatrixArray, this->MatrixArray, 4, 4);
    
    return ResultMatrix;
}

SPTH::A3DMatrix SPTH::A4DMatrix::Cut (int Row, int Column) const
{
    A3DMatrix ResultMatrix;
    int Index = 0;
    
    for (int RowIndex = 0; RowIndex < 4; ++RowIndex)
    {
        for (int ColumnIndex = 0; ColumnIndex < 4; ++ColumnIndex)
        {
            if (RowIndex == Row || ColumnIndex == Column)
            {
                continue;
            }
            ResultMatrix.MatrixArray[Index++] = this->MatrixArray[4 * RowIndex + ColumnIndex];
        }
    }
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::A4DMatrix::Minor () const
{
    A4DMatrix ResultMatrix;
    
    for (int RowIndex = 0; RowIndex < 4; ++RowIndex)
    {
        for (int ColumnIndex = 0; ColumnIndex < 4; ++ColumnIndex)
        {
            ResultMatrix[RowIndex][ColumnIndex] = this->Cut(RowIndex, ColumnIndex).Determinant();
        }
    }
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::A4DMatrix::Cofactor () const
{
    A4DMatrix ResultMatrix;
    SPTH::Cofactor(ResultMatrix.MatrixArray, this->Minor().MatrixArray, 4, 4);
    
    return ResultMatrix;
}

float SPTH::A4DMatrix::Determinant () const
{
    float ResultDeterminant = 0.0f;
    A4DMatrix ResultCofactorMatrix = this->Cofactor();
    
    for (int Index = 0; Index < 4; ++Index)
    {
        ResultDeterminant += this->MatrixArray[Index] * ResultCofactorMatrix[0][Index];
    }
    
    return ResultDeterminant;
}

SPTH::A4DMatrix SPTH::A4DMatrix::Adjugate () const
{
    return this->Cofactor().Transpose();
}

SPTH::A4DMatrix SPTH::A4DMatrix::Inverse () const
{
    float ResultDeterminant = this->Determinant();
    
    if (COMPARISON(ResultDeterminant, 0.0f))
    {
        return *this;
    }
    
    return this->Adjugate() * (1.0f / ResultDeterminant);
}

SPTH::A3DVector SPTH::A4DMatrix::GetTranslation () const
{
    return {this->M41, this->M42, this->M43};
}

SPTH::A3DVector SPTH::A4DMatrix::GetScale () const
{
    return {this->M11, this->M22, this->M33};
}

SPTH::A3DVector SPTH::A4DMatrix::MultiplyPoint (const A3DVector &Point) const
{
    A3DVector ResultVector;
    ResultVector.X = Point.X * this->M11 + Point.Y * this->M21 + Point.Z * this->M31 + this->M41;
    ResultVector.Y = Point.X * this->M12 + Point.Y * this->M22 + Point.Z * this->M32 + this->M42;
    ResultVector.Z = Point.X * this->M13 + Point.Y * this->M23 + Point.Z * this->M33 + this->M43;
    
    return ResultVector;
}

SPTH::A3DVector SPTH::A4DMatrix::MultiplyVector (const A3DVector &Vector) const
{
    A3DVector ResultVector;
    ResultVector.X = Vector.X * this->M11 + Vector.Y * this->M21 + Vector.Z * this->M31;
    ResultVector.Y = Vector.X * this->M12 + Vector.Y * this->M22 + Vector.Z * this->M32;
    ResultVector.Z = Vector.X * this->M13 + Vector.Y * this->M23 + Vector.Z * this->M33;
    
    return ResultVector;
}

void SPTH::Transpose (float *ResultMatrixArray, const float *Matrix, int Row, int Column)
{
    for (int Index = 0; Index < Row * Column; ++Index) {
        int Row_ = Index / Row;
        int Column_ = Index % Row;
        ResultMatrixArray[Index] = Matrix[Row * Column_ + Row_];
    }
}

void SPTH::Multiply (float *ResultMatrixArray, const float *MatrixA, int MatrixARow, int MatrixAColumn, const float *MatrixB, int MatrixBRow, int MatrixBColumn)
{
    assert (MatrixAColumn == MatrixBRow);
    for (int i = 0; i < MatrixARow; ++i)
    {
        for (int j = 0; j < MatrixBColumn; ++j)
        {
            ResultMatrixArray[MatrixBColumn * i + j] = 0.0f;
            for (int k = 0; k < MatrixBRow; ++k)
            {
                int IndexMatrixA = MatrixAColumn * i + k;
                int IndexMatrixB = MatrixBColumn * k + j;
                ResultMatrixArray[MatrixBColumn * i + j] += MatrixA[IndexMatrixA] * MatrixB[IndexMatrixB];
            }
        }
    }
    
    /*
    unsigned int ResultIndex = 0, Total = MatrixARows * MatrixAColumns;
    for (int i = 0; i < Total; i += MatrixAColumns) {
        for (int j = 0; j < MatrixBRows; ++j) {
            int A = j;
            ResultMatrix[ResultIndex] = 0.0f;
            for (int k = 0; k < MatrixARows; ++k) {
                ResultMatrix[ResultIndex] += MatrixA[k + i] * MatrixB[(k / MatrixARows) + (A % (Total))];
                A += MatrixBRows;
            }
            ResultIndex++;
        }
    }
    */
}

void SPTH::Cofactor (float *ResultMatrixArray, const float *Matrix, int Row, int Column)
{
    for (int RowIndex = 0; RowIndex < Row; ++RowIndex)
    {
        for (int ColumnIndex = 0; ColumnIndex < Column; ++ColumnIndex)
        {
            int To = Column * ColumnIndex + RowIndex;
            int From = Column * ColumnIndex + RowIndex;
            float Sign = powf(-1.0f, (float)(RowIndex + ColumnIndex));
            ResultMatrixArray[To] = Sign * Matrix[From];
        }
    }
}

SPTH::A3DMatrix SPTH::Rotation3D (float XAngle, float YAngle, float ZAngle)
{
    return ZRotation3D(ZAngle) * XRotation3D(XAngle) * YRotation3D(YAngle);
}

SPTH::A3DMatrix SPTH::ZRotation3D (float ZAngle)
{
    ZAngle = DEGREES_TO_RADIANS(ZAngle);
    
    return {cosf(ZAngle), sinf(ZAngle), 0.0f,
            -sinf(ZAngle), cosf(ZAngle), 0.0f,
            0.0f, 0.0f, 1.0f};
}

SPTH::A3DMatrix SPTH::XRotation3D (float XAngle)
{
    XAngle = DEGREES_TO_RADIANS(XAngle);
    
    return {1.0f, 0.0f, 0.0f,
            0.0f, cosf(XAngle), sinf(XAngle),
            0.0f, -sinf(XAngle), cos(XAngle)};
}

SPTH::A3DMatrix SPTH::YRotation3D (float YAngle)
{
    YAngle = DEGREES_TO_RADIANS(YAngle);
    
    return {cosf(YAngle), 0.0f, -sinf(YAngle),
            0.0f, 1.0f, 0.0f,
            sinf(YAngle), 0.0f, cosf(YAngle)};
}

SPTH::A3DMatrix SPTH::AxisAngle3D (const A3DVector &Axis, float Angle)
{
    Angle = DEGREES_TO_RADIANS(Angle);
    float Cosine = cosf(Angle);
    float Sine = sinf(Angle);
    float T = 1.0f - cosf(Angle);
    float X = Axis.X;
    float Y = Axis.Y;
    float Z = Axis.Z;
    
    if (!COMPARISON(Axis.MagnitudeSquared(), 1.0f))
    {
        float InverseMagnitude = 1.0f / Axis.Magnitude();
        X *= InverseMagnitude;
        Y *= InverseMagnitude;
        Z *= InverseMagnitude;
    }
    
    return {T * (X * X) + Cosine, T * X * Y + Sine * Z, T * X * Z + Sine * Y,
            T * X * Y - Sine * Z, T * (Y * Y) + Cosine, T * Y * Z + Sine * X,
            T * X * Z + Sine * Y, T * Y * Z - Sine * X, T * (Z * Z) + Cosine};
}

SPTH::A4DMatrix SPTH::Translation (float X, float Y, float Z)
{
    A4DMatrix ResultMatrix;
    ResultMatrix.M41 = X;
    ResultMatrix.M42 = Y;
    ResultMatrix.M43 = Z;
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::Translation (const A3DVector &TranslationVector)
{
    A4DMatrix ResultMatrix;
    ResultMatrix.M41 = TranslationVector.X;
    ResultMatrix.M42 = TranslationVector.Y;
    ResultMatrix.M43 = TranslationVector.Z;
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::FromMatrix3D (const A3DMatrix &Matrix)
{
    A4DMatrix Result;
    Result.M11 = Matrix.M11;
    Result.M12 = Matrix.M12;
    Result.M13 = Matrix.M13;
    Result.M21 = Matrix.M21;
    Result.M22 = Matrix.M22;
    Result.M23 = Matrix.M23;
    Result.M31 = Matrix.M31;
    Result.M32 = Matrix.M32;
    Result.M33 = Matrix.M33;
    
    return Result;
}

SPTH::A4DMatrix SPTH::Scale (float X, float Y, float Z)
{
    A4DMatrix ResultMatrix;
    ResultMatrix.M11 = X;
    ResultMatrix.M22 = Y;
    ResultMatrix.M33 = Z;
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::Scale (const A3DVector &ScaleVector)
{
    A4DMatrix ResultMatrix;
    ResultMatrix.M11 = ScaleVector.X;
    ResultMatrix.M22 = ScaleVector.Y;
    ResultMatrix.M33 = ScaleVector.Z;
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::Rotation4D (float XAngle, float YAngle, float ZAngle)
{
    return ZRotation4D(ZAngle) * XRotation4D(XAngle) * YRotation4D(YAngle);
}

SPTH::A4DMatrix SPTH::ZRotation4D (float ZAngle)
{
    ZAngle = DEGREES_TO_RADIANS(ZAngle);
    
    return {cosf(ZAngle), sinf(ZAngle), 0.0f, 0.0f,
            -sinf(ZAngle), cosf(ZAngle), 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f};
}

SPTH::A4DMatrix SPTH::XRotation4D (float XAngle)
{
    XAngle = DEGREES_TO_RADIANS(XAngle);
    
    return {1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, cosf(XAngle), sinf(XAngle), 0.0f,
            0.0f, -sinf(XAngle), cos(XAngle), 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f};
}

SPTH::A4DMatrix SPTH::YRotation4D (float YAngle)
{
    YAngle = DEGREES_TO_RADIANS(YAngle);
    
    return {cosf(YAngle), 0.0f, -sinf(YAngle), 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            sinf(YAngle), 0.0f, cosf(YAngle), 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f};
}

SPTH::A4DMatrix SPTH::AxisAngle4D (const A3DVector &Axis, float Angle)
{
    Angle = DEGREES_TO_RADIANS(Angle);
    float Cosine = cosf(Angle);
    float Sine = sinf(Angle);
    float T = 1.0f - cosf(Angle);
    float X = Axis.X;
    float Y = Axis.Y;
    float Z = Axis.Z;
    
    if (!COMPARISON(Axis.MagnitudeSquared(), 1.0f))
    {
        float InverseMagnitude = 1.0f / Axis.Magnitude();
        X *= InverseMagnitude;
        Y *= InverseMagnitude;
        Z *= InverseMagnitude;
    }
    
    return {T * (X * X) + Cosine, T * X * Y + Sine * Z, T * X * Z + Sine * Y, 0.0f,
            T * X * Y - Sine * Z, T * (Y * Y) + Cosine, T * Y * Z + Sine * X, 0.0f,
            T * X * Z + Sine * Y, T * Y * Z - Sine * X, T * (Z * Z) + Cosine, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f};
}

SPTH::A4DMatrix SPTH::Transformation (const A3DVector &ScaleVector, const A3DVector &EulerRotationVector, const A3DVector &TranslationVector)
{
    return Scale(ScaleVector) * Rotation4D(EulerRotationVector.X, EulerRotationVector.Y, EulerRotationVector.Z) * Translation(TranslationVector);
}

SPTH::A4DMatrix SPTH::Transformation (const A3DVector &ScaleVector, const A3DVector &AxisRotationVector, float Angle, const A3DVector &TranslationVector)
{
    return Scale(ScaleVector) * AxisAngle4D(AxisRotationVector, Angle) * Translation(TranslationVector);
}

SPTH::A4DMatrix SPTH::LookAt (const A3DVector &CameraPosition, const A3DVector &TargetPosition, const A3DVector &UpwardVector)
{
    A3DVector ForwardVector = (TargetPosition - CameraPosition).Normalized();
    A3DVector RightVector = (ForwardVector.Cross(UpwardVector)).Normalized();
    A3DVector EnsureUpVector = ForwardVector.Cross(RightVector);
    
    return {RightVector.X, EnsureUpVector.X, ForwardVector.X, 0.0f,
            RightVector.Y, EnsureUpVector.Y, ForwardVector.Y, 0.0f,
            RightVector.Z, EnsureUpVector.Z,ForwardVector.Z, 0.0f,
            -RightVector.Dot(CameraPosition),
            -EnsureUpVector.Dot(CameraPosition),
            -ForwardVector.Dot(CameraPosition), 1.0f};
}

SPTH::A4DMatrix SPTH::Perspective (float FOV, float Aspect, float ZNear, float ZFar)
{
    float TangentOfHalfFOV = tanf(DEGREES_TO_RADIANS(FOV * 0.5f));
    float Cotangent = 1.0f / TangentOfHalfFOV;
    A4DMatrix ResultMatrix;
    ResultMatrix.M11 = Cotangent / Aspect;
    ResultMatrix.M22 = Cotangent;
    ResultMatrix.M33 = ZFar / (ZFar - ZNear);
    ResultMatrix.M34 = 1.0f;
    ResultMatrix.M43 = -ZNear * ResultMatrix.M33;
    ResultMatrix.M44 = 0.0f;
    
    return ResultMatrix;
}

SPTH::A4DMatrix SPTH::Orthographic (float Left, float Right, float Bottom, float Top, float ZNear, float ZFar)
{
    float M11 = 2.0f / (Right - Left);
    float M22 = 2.0f / (Top - Bottom);
    float M33 = 1.0f / (ZFar - ZNear);
    float M41 = (Left * Right) / (Left - Right);
    float M42 = (Top + Bottom) / (Bottom - Top);
    float M43 = ZNear / (ZNear - ZFar);
    
    return {M11, 0.0f, 0.0f, 0.0f,
            0.0f, M22, 0.0f, 0.0f,
            0.0f, 0.0f, M33, 0.0f,
            M41, M42, M43, 1.0f};
}

#ifndef NO_EXTRAS
std::ostream& operator << (std::ostream &OS, const SPTH::A2DMatrix &MatrixA)
{
    OS << "[" << MatrixA.M11 << "\t" << MatrixA.M12 << std::endl;
    OS << MatrixA.M21 << "\t" << MatrixA.M22 << "]";
    
    return OS;
}

std::ostream& operator << (std::ostream &OS, const SPTH::A3DMatrix &MatrixA)
{
    OS << "[" << MatrixA.M11 << "\t" << MatrixA.M12 << "\t" << MatrixA.M13 << std::endl;
    OS << MatrixA.M21 << "\t" << MatrixA.M22 << "\t" << MatrixA.M23 << std::endl;
    OS << MatrixA.M31 << "\t" << MatrixA.M32 << "\t" << MatrixA.M33 << "]";
    
    return OS;
}

std::ostream& operator << (std::ostream &OS, const SPTH::A4DMatrix &MatrixA)
{
    OS << "[" << MatrixA.M11 << "\t" << MatrixA.M12 << "\t" << MatrixA.M13 << "\t" << MatrixA.M14 << std::endl;
    OS << MatrixA.M21 << "\t" << MatrixA.M22 << "\t" << MatrixA.M23 << "\t" << MatrixA.M24 << std::endl;
    OS << MatrixA.M31 << "\t" << MatrixA.M32 << "\t" << MatrixA.M33 << "\t" << MatrixA.M34 << std::endl;
    OS << MatrixA.M41 << "\t" << MatrixA.M42 << "\t" << MatrixA.M43 << "\t" << MatrixA.M44 << "]";
    
    return OS;
}
#endif