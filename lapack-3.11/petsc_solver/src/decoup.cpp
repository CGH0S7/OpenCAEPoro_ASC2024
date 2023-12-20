#include "PETScSolver.h"

#include "cblas.h"
#include "lapacke.h"

void inverse(double *A, int N)
{
    int *IPIV = new int[N];
    int LWORK = N * N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N, &N, A, &N, IPIV, &INFO);

    if (INFO > 0)
    {
        printf("WARNING: The factorization has been completed, but the factor U is exactly singular");
    }
    else if (INFO < 0)
    {
        printf("WARNING: There is an illegal value");
    }

    dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

    delete[] IPIV;
    delete[] WORK;
}

void smat_mul(double *A, double *B, double *C, int nb)
{
    const int n = nb;
    // 全0 初始化 c
    memset(C, 0, n * n * sizeof(double));
    int i, j, k;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            for (k = 0; k < n; ++k)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void smat_identity(double *A, int nb)
{
    const int n = nb;
    // 全0 初始化 A
    memset(A, 0, n * n * sizeof(double));
    int i;
    for (i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }
}

void smat_vec_mul(double *A, double *b, double *c, int nb)
{
    const int n = nb;
    // 全0 初始化 c
    memset(c, 0, n * sizeof(double));
    int i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            c[i] += A[i * n + j] * b[j];
}

void smat_identity_4x4(double *A)
{
    const int n = 4;
    // 全0 初始化 A
    memset(A, 0, n * n * sizeof(double));
    int i;
    for (i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }
}

void smat_identity_2x2(double *A)
{
    const int n = 2;
    // 全0 初始化 A
    memset(A, 0, n * n * sizeof(double));
    int i;
    for (i = 0; i < n; ++i)
    {
        A[i * n + i] = 1.0;
    }
}

void smat_mul_4x4(double *A, double *B, double *C)
{
    const int n = 4;
    // 全0 初始化 c
    memset(C, 0, n * n * sizeof(double));
    int i, j, k;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            for (k = 0; k < n; ++k)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void smat_mul_2x2(double *A, double *B, double *C)
{
    const int n = 2;
    // 全0 初始化 c
    memset(C, 0, n * n * sizeof(double));
    int i, j, k;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            for (k = 0; k < n; ++k)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void smat_vec_mul_4(double *A, double *b, double *c)
{
    const int n = 4;
    // 全0 初始化 c
    memset(c, 0, n * sizeof(double));
    int i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            c[i] += A[i * n + j] * b[j];
}

void smat_vec_mul_2(double *A, double *b, double *c)
{
    const int n = 2;
    // 全0 初始化 c
    memset(c, 0, n * sizeof(double));
    int i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            c[i] += A[i * n + j] * b[j];
}

void smat_inv_4x4(double *a)
{
    const double a11 = a[0], a12 = a[1], a13 = a[2], a14 = a[3];
    const double a21 = a[4], a22 = a[5], a23 = a[6], a24 = a[7];
    const double a31 = a[8], a32 = a[9], a33 = a[10], a34 = a[11];
    const double a41 = a[12], a42 = a[13], a43 = a[14], a44 = a[15];

    const double M11 = a22 * a33 * a44 + a23 * a34 * a42 + a24 * a32 * a43 - a22 * a34 * a43 - a23 * a32 * a44 - a24 * a33 * a42;
    const double M12 = a12 * a34 * a43 + a13 * a32 * a44 + a14 * a33 * a42 - a12 * a33 * a44 - a13 * a34 * a42 - a14 * a32 * a43;
    const double M13 = a12 * a23 * a44 + a13 * a24 * a42 + a14 * a22 * a43 - a12 * a24 * a43 - a13 * a22 * a44 - a14 * a23 * a42;
    const double M14 = a12 * a24 * a33 + a13 * a22 * a34 + a14 * a23 * a32 - a12 * a23 * a34 - a13 * a24 * a32 - a14 * a22 * a33;
    const double M21 = a21 * a34 * a43 + a23 * a31 * a44 + a24 * a33 * a41 - a21 * a33 * a44 - a23 * a34 * a41 - a24 * a31 * a43;
    const double M22 = a11 * a33 * a44 + a13 * a34 * a41 + a14 * a31 * a43 - a11 * a34 * a43 - a13 * a31 * a44 - a14 * a33 * a41;
    const double M23 = a11 * a24 * a43 + a13 * a21 * a44 + a14 * a23 * a41 - a11 * a23 * a44 - a13 * a24 * a41 - a14 * a21 * a43;
    const double M24 = a11 * a23 * a34 + a13 * a24 * a31 + a14 * a21 * a33 - a11 * a24 * a33 - a13 * a21 * a34 - a14 * a23 * a31;
    const double M31 = a21 * a32 * a44 + a22 * a34 * a41 + a24 * a31 * a42 - a21 * a34 * a42 - a22 * a31 * a44 - a24 * a32 * a41;
    const double M32 = a11 * a34 * a42 + a12 * a31 * a44 + a14 * a32 * a41 - a11 * a32 * a44 - a12 * a34 * a41 - a14 * a31 * a42;
    const double M33 = a11 * a22 * a44 + a12 * a24 * a41 + a14 * a21 * a42 - a11 * a24 * a42 - a12 * a21 * a44 - a14 * a22 * a41;
    const double M34 = a11 * a24 * a32 + a12 * a21 * a34 + a14 * a22 * a31 - a11 * a22 * a34 - a12 * a24 * a31 - a14 * a21 * a32;
    const double M41 = a21 * a33 * a42 + a22 * a31 * a43 + a23 * a32 * a41 - a21 * a32 * a43 - a22 * a33 * a41 - a23 * a31 * a42;
    const double M42 = a11 * a32 * a43 + a12 * a33 * a41 + a13 * a31 * a42 - a11 * a33 * a42 - a12 * a31 * a43 - a13 * a32 * a41;
    const double M43 = a11 * a23 * a42 + a12 * a21 * a43 + a13 * a22 * a41 - a11 * a22 * a43 - a12 * a23 * a41 - a13 * a21 * a42;
    const double M44 = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31;

    const double det = a11 * M11 + a12 * M21 + a13 * M31 + a14 * M41;
    double det_inv;

    if (fabs(det) < 1e-22)
    {
        printf("### WARNING: 4x4 Matrix is nearly singular! det = %e\n", det);
    }

    det_inv = 1.0 / det;

    a[0] = M11 * det_inv;
    a[1] = M12 * det_inv;
    a[2] = M13 * det_inv;
    a[3] = M14 * det_inv;
    a[4] = M21 * det_inv;
    a[5] = M22 * det_inv;
    a[6] = M23 * det_inv;
    a[7] = M24 * det_inv;
    a[8] = M31 * det_inv;
    a[9] = M32 * det_inv;
    a[10] = M33 * det_inv;
    a[11] = M34 * det_inv;
    a[12] = M41 * det_inv;
    a[13] = M42 * det_inv;
    a[14] = M43 * det_inv;
    a[15] = M44 * det_inv;
}

void smat_inv_2x2(double *a)
{
    const double a11 = a[0], a12 = a[1], a21 = a[2], a22 = a[3];

    const double M11 = a22;
    const double M12 = -a12;
    const double M21 = -a21;
    const double M22 = a11;

    const double det = a11 * a22 - a12 * a21;
    double det_inv;

    if (fabs(det) < 1e-22)
    {
        printf("### WARNING: 2x2 Matrix is nearly singular! det = %e\n", det);
    }

    det_inv = 1.0 / det;

    a[0] = M11 * det_inv;
    a[1] = M12 * det_inv;
    a[2] = M21 * det_inv;
    a[3] = M22 * det_inv;
}

// ABF
void decoup_abf_2x2(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow)
{
    int nb2 = nb * nb;
    double *smat = (double *)malloc(nb2 * sizeof(double));
    double *binv = (double *)malloc(nb2 * sizeof(double));
    double *svec = (double *)malloc(nb * sizeof(double));

    int i, j, ibegin, iend;
    for (i = 0; i < nrow; ++i)
    {
        // 将 val 拷贝给 binv
        memcpy(binv, val, nb2 * sizeof(double));
        // 求 binv 的逆
        smat_inv_2x2(binv);
        ibegin = rpt[i];
        iend = rpt[i + 1];
        for (j = ibegin; j < iend; ++j)
        {
            if (j != ibegin)
            {
                // 矩阵乘 smat = binv * val
                smat_mul_2x2(binv, val, smat);
                // 将 smat 拷贝给 val
                memcpy(val, smat, nb2 * sizeof(double));
            }
            else
            {
                // 矩阵单位化
                smat_identity_2x2(val);
            }
            val += nb2;
        }
        // 矩阵乘向量 smat = binv * rhs
        smat_vec_mul_2(binv, rhs, smat);
        // 将 smat 拷贝给 rhs
        memcpy(rhs, smat, nb * sizeof(double));
        rhs += nb;
    }

    free(smat);
    free(binv);
    free(svec);
}

// ABF
void decoup_abf_4x4(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow)
{
    int nb2 = nb * nb;
    double *smat = (double *)malloc(nb2 * sizeof(double));
    double *binv = (double *)malloc(nb2 * sizeof(double));
    double *svec = (double *)malloc(nb * sizeof(double));

    int i, j, ibegin, iend;
    for (i = 0; i < nrow; ++i)
    {
        // 将 val 拷贝给 binv
        memcpy(binv, val, nb2 * sizeof(double));
        // 求 binv 的逆
        smat_inv_4x4(binv);
        ibegin = rpt[i];
        iend = rpt[i + 1];
        for (j = ibegin; j < iend; ++j)
        {
            if (j != ibegin)
            {
                // 矩阵乘 smat = binv * val
                smat_mul_4x4(binv, val, smat);
                // 将 smat 拷贝给 val
                memcpy(val, smat, nb2 * sizeof(double));
            }
            else
            {
                // 矩阵单位化
                smat_identity_4x4(val);
            }
            val += nb2;
        }
        // 矩阵乘向量 smat = binv * rhs
        smat_vec_mul_4(binv, rhs, smat);
        // 将 smat 拷贝给 rhs
        memcpy(rhs, smat, nb * sizeof(double));
        rhs += nb;
    }

    free(smat);
    free(binv);
    free(svec);
}

// ABF
void decoup_abf_nb(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow)
{
    int nb2 = nb * nb;
    double *smat = (double *)malloc(nb2 * sizeof(double));
    double *binv = (double *)malloc(nb2 * sizeof(double));
    double *svec = (double *)malloc(nb * sizeof(double));

    int i, j, ibegin, iend;
    for (i = 0; i < nrow; ++i)
    {
        // 将 val 拷贝给 binv
        memcpy(binv, val, nb2 * sizeof(double));
        // 求 binv 的逆
        inverse(binv, nb);
        ibegin = rpt[i];
        iend = rpt[i + 1];
        for (j = ibegin; j < iend; ++j)
        {
            if (j != ibegin)
            {
                // 矩阵乘 smat = binv * val
                smat_mul(binv, val, smat, nb);
                // 将 smat 拷贝给 val
                memcpy(val, smat, nb2 * sizeof(double));
            }
            else
            {
                // 矩阵单位化
                smat_identity(val, nb);
            }
            val += nb2;
        }
        // 矩阵乘向量 smat = binv * rhs
        smat_vec_mul(binv, rhs, smat, nb);
        // 将 smat 拷贝给 rhs
        memcpy(rhs, smat, nb * sizeof(double));
        rhs += nb;
    }

    free(smat);
    free(binv);
    free(svec);
}

// ABF decoupling method
void decouple_abf(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow)
{
    switch (nb)
    {
    case 2:
        decoup_abf_2x2(val, rhs, rpt, cpt, nb, nrow);
        break;

    case 4:
        decoup_abf_4x4(val, rhs, rpt, cpt, nb, nrow);
        break;

    default:
        decoup_abf_nb(val, rhs, rpt, cpt, nb, nrow);
        break;
    }
}

// Quasi-IMPES decoupling method
void decouple_QI(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow, int is_thermal)
{
    int nb2 = nb * nb;
    double *smat = (double *)malloc(nb2 * sizeof(double));
    double *binv = (double *)malloc(nb2 * sizeof(double));
    double *svec = (double *)malloc(nb * sizeof(double));

    int i, j, l, ibegin, iend;
    for (i = 0; i < nrow; ++i)
    {
        // form the diagonal sub-blocks for analytical decoupling
        smat_identity(binv, nb);
        if (is_thermal)
        {
            /**
             * A =
             * [PP  PN1  PN2  ... PT ]
             * [N1P N1N1 N1N2 ... N1T]
             * [N2P N2N1 N2N2 ... N2T]
             *         ...
             * [TP  TN1  TN2  ... TT]
             * */
            for (l = 0; l < nb - 2; l++)
            {
                binv[1 + l] = -val[1 + l] / val[(l + 1) * nb + l + 1];
                binv[(nb - 1) * nb + 1 + l] = -val[(nb - 1) * nb + 1 + l] / val[(l + 1) * nb + l + 1];
            }
        }
        else
        {
            for (l = 0; l < nb - 1; l++)
                binv[1 + l] = -val[1 + l] / val[(l + 1) * nb + l + 1];
        }

        // compute D^{-1}*A
        ibegin = rpt[i];
        iend = rpt[i + 1];
        for (j = ibegin; j < iend; ++j)
        {
            // 矩阵乘 smat = binv * val
            smat_mul(binv, val, smat, nb);
            // 将 smat 拷贝给 val
            memcpy(val, smat, nb2 * sizeof(double));

            val += nb2;
        }
        // 矩阵乘向量 smat = binv * rhs
        smat_vec_mul(binv, rhs, smat, nb);
        // 将 smat 拷贝给 rhs
        memcpy(rhs, smat, nb * sizeof(double));
        rhs += nb;
    }

    free(smat);
    free(binv);
    free(svec);
}

// Analytical decoupling method
void decouple_anl(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow, int is_thermal)
{
    int nb2 = nb * nb;
    double *smat = (double *)malloc(nb2 * sizeof(double));
    double *binv = (double *)malloc(nb2 * sizeof(double));
    double *svec = (double *)malloc(nb * sizeof(double));

    int i, j, l, ibegin, iend;
    for (i = 0; i < nrow; ++i)
    {
        // form the diagonal sub-blocks for analytical decoupling
        smat_identity(binv, nb);
        if (is_thermal)
        {
            /**
             * A =
             * [PP  PN1  PN2  ... PT ]
             * [N1P N1N1 N1N2 ... N1T]
             * [N2P N2N1 N2N2 ... N2T]
             *         ...
             * [TP  TN1  TN2  ... TT]
             * */
            for (l = 0; l < nb - 2; l++)
            {
                binv[1 + l] = -val[1 + l];
                // binv[nb2 - 2 - l] = -val[nb2 - 2 - l];
                binv[(nb - 1) * nb + 1 + l] = -val[(nb - 1) * nb + 1 + l];
            }
        }
        else
        {
            for (l = 0; l < nb - 1; l++)
                binv[1 + l] = -val[1 + l];
        }

        // compute D^{-1}*A
        ibegin = rpt[i];
        iend = rpt[i + 1];
        for (j = ibegin; j < iend; ++j)
        {
            // 矩阵乘 smat = binv * val
            smat_mul(binv, val, smat, nb);
            // 将 smat 拷贝给 val
            memcpy(val, smat, nb2 * sizeof(double));

            val += nb2;
        }
        // 矩阵乘向量 smat = binv * rhs
        smat_vec_mul(binv, rhs, smat, nb);
        // 将 smat 拷贝给 rhs
        memcpy(rhs, smat, nb * sizeof(double));
        rhs += nb;
    }

    free(smat);
    free(binv);
    free(svec);
}

// Semi-analytical decoupling method
void decouple_sem(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow, int is_thermal)
{
    int nb2 = nb * nb;
    double *smat = (double *)malloc(nb2 * sizeof(double));
    double *binv = (double *)malloc(nb2 * sizeof(double));
    double *svec = (double *)malloc(nb * sizeof(double));

    int i, j, l, ibegin, iend;
    for (i = 0; i < nrow; ++i)
    {
        // S1: Form the ABF part first
        memcpy(binv, val, nb2 * sizeof(double));
        inverse(binv, nb);

        // S2: Replace the first line with analytical decoupling
        if (is_thermal)
        {
            /**
             * A =
             * [PP  PN1  PN2  ... PT ]
             * [N1P N1N1 N1N2 ... N1T]
             * [N2P N2N1 N2N2 ... N2T]
             *         ...
             * [TP  TN1  TN2  ... TT]
             * */
            // Replace the first and final lines with analytical decoupling
            binv[0] = 1;
            binv[nb - 1] = 0;
            binv[(nb - 1) * nb] = 0;
            binv[nb2 - 1] = 1;
            for (l = 0; l < nb - 2; l++)
            {
                binv[1 + l] = -val[1 + l];
                // binv[nb2 - 2 - l] = -val[nb2 - 2 - l];
                binv[(nb - 1) * nb + 1 + l] = -val[(nb - 1) * nb + 1 + l];
            }
        }
        else
        {
            // Replace the first line with analytical decoupling
            binv[0] = 1;
            for (l = 0; l < nb - 1; l++)
                binv[1 + l] = -val[1 + l];
        }

        // compute D^{-1}*A
        ibegin = rpt[i];
        iend = rpt[i + 1];
        for (j = ibegin; j < iend; ++j)
        {
            // 矩阵乘 smat = binv * val
            smat_mul(binv, val, smat, nb);
            // 将 smat 拷贝给 val
            memcpy(val, smat, nb2 * sizeof(double));

            val += nb2;
        }
        // 矩阵乘向量 smat = binv * rhs
        smat_vec_mul(binv, rhs, smat, nb);
        // 将 smat 拷贝给 rhs
        memcpy(rhs, smat, nb * sizeof(double));
        rhs += nb;
    }

    free(smat);
    free(binv);
    free(svec);
}

void decoup(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow, int decoup_type, int is_thermal)
{
    switch (decoup_type)
    {
    case 0: // don't using decoupling method (None)
            // do nothing
        break;

    case 2: // Analytical decoupling method (ANL)
        decouple_anl(val, rhs, rpt, cpt, nb, nrow, is_thermal);
        break;

    case 3: // Semi-analytical decoupling method (SEM)
        decouple_sem(val, rhs, rpt, cpt, nb, nrow, is_thermal);
        break;

    case 4: // Quasi-IMPES decoupling method (QI)
        decouple_QI(val, rhs, rpt, cpt, nb, nrow, is_thermal);
        break;

    default: // case 1: ABF
        decouple_abf(val, rhs, rpt, cpt, nb, nrow);
        break;
    }
}
