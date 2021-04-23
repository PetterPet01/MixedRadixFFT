using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Lomont;

namespace MixedRadixFFTTest1
{
    public class MixedFFT
    {
        LomontFFT fft;
        int N;
        double[][] splitData;
        double[,] twF;
        int M, P;
        public MixedFFT(int M, int P, int A = 1, int B = -1)
        {
            fft = new LomontFFT();
            fft.A = A;
            fft.B = B;
            this.M = M;
            this.P = P;
            N = M * P;
            int i, j;
            splitData = new double[M][];
            
            twF = new double[N + 2, M];
            for (i = 0; i < N + 2; i += 2)
                for (j = 0; j < M; j++)
                {
                    twF[i, j] = Math.Cos(Math.PI * i * j / N);
                    twF[i + 1, j] = -Math.Sin(Math.PI * i * j / N);
                }
        }
        public double[] Compute(double[] realData)
        {
            if (M * P != realData.Length)
                throw new Exception();
            int i, j;
            for (i = 0; i < M; i++)
                splitData[i] = new double[P];
            for (i = 0; i < M; i++)
                for (j = 0; j < P; j++)
                    splitData[i][j] = realData[i + j * M];
            for (i = 0; i < M; i++)
            {
                fft.RealFFT(splitData[i], true);
                splitData[i] = convertToCmplxFFT(splitData[i]);
            }
            int irsx1 = 0; //1st special real-valued only index for X (fft-computed split data)
            int irsx2 = N % (P * 2); //2nd special real-valued only index for X (fft-computed split data)
            double[] output = new double[N];
            for (j = 1; j < M; j++)
            {
                for (i = 2; i < N; i += 2)
                {
                    int ir = i;
                    int ii = i + 1;
                    int irx = ir % (P * 2);
                    int iix = ii % (P * 2);

                    output[ir] += splitData[j][irx] * twF[ir, j] - splitData[j][iix] * twF[ii, j];
                    output[ii] += splitData[j][irx] * twF[ii, j] + splitData[j][iix] * twF[ir, j];
                }
                output[0] += splitData[j][irsx1] * twF[0, j];
                output[1] += splitData[j][irsx2] * twF[N, j];
            }
            for (i = 2; i < N; i += 2)
            {
                int ir = i;
                int ii = i + 1;
                int irx = ir % (P * 2);
                int iix = ii % (P * 2);

                output[ir] += splitData[0][irx];
                output[ii] += splitData[0][iix];
            }
            output[0] += splitData[0][irsx1] * twF[0, 0];
            output[1] += splitData[0][irsx2] * twF[N, 0];
            return output;
        }

        double[] convertToCmplxFFT(double[] data)
        {
            int datLen = data.Length;
            double[] output = new double[datLen * 2];
            output[0] = data[0];
            output[datLen] = data[1];
            int uniqueDatLen = datLen - 2;
            Array.Copy(data, 2, output, 2, uniqueDatLen);
            int i;
            int j = uniqueDatLen;
            for (i = datLen + 2; i < datLen * 2; i += 2, j -= 2)
            {
                output[i] = data[j];
                output[i + 1] = -data[j + 1];
            }
            return output;
        }
    }
}
