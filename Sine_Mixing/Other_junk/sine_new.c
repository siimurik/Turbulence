// if N is a power of 2 then N = 0100000b=32
// then N-1 = 0011111b
// then i = (i + int(a*np.sin(j*2*np.pi/N))) & N-1

// 0011111 & 
// 0111101 = 29+32
// 0011101 = 29          bitwise and
// i = (i + a*np.sin(j*2*np.pi/N))
// 000000b-1 = 111111b
// conclusion: n % 32 = n & 31
#define N 100

int main(void){
    int i,j;
    static double c_new[N][N], c[N][N];
    float a;

    a = 1.0;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            c_new[i + (int)(a*sin(j*2*pi/N)) & N-1][j] = c[i][j];
        }
    }
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            c_new[j][i + (int)(a*sin(j*2*pi/N)) & N-1] = c[i][j];
        }
    }

  return 0;
}