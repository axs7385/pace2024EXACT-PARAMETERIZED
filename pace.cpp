#include <bits/stdc++.h>
#include <signal.h>
#include <unistd.h>

#define sig

#define exact

using namespace std;
#define ll long long
#ifdef sig
volatile sig_atomic_t tle = 0;
void term(int signum)
{
    tle = 1;
}
#else
int tle = 0;
#endif
const int mod = 998244353, base = 19260817; // hash 常数

set<ll> tabu;                    // 禁忌表
const int N = 300005, M = 21000; // 总店
int n, m, K;
vector<int> e[N];
// int c[M][M],p[M+1][M+1];
int ord[N], _ord[N];
int ords[10][M + 5], used[M + 5], child[M + 5];
ll ans[10], childans;
vector<ll> c[M], p[M + 1];
double avg[N];

// 升序排列
bool cmp(int x, int y)
{
    return avg[x] < avg[y];
}

// 读入数据
void read(int &x)
{
    x = 0;
    char c = getchar();
    while (c == 'c')
    {
        while (c != '\n')
            c = getchar();
        c = getchar();
    }
    while (c > '9' || c < '0')
        c = getchar();
    while (c >= '0' && c <= '9')
    {
        x = x * 10 + c - '0';
        c = getchar();
    }
}

ll calc(int i, int j)
{
    int l = 0;
    int lj = e[j].size();
    ll sum = 0;
    for (int x : e[i])
    {
        while (l < lj && e[j][l] < x)
            ++l;
        sum += l;
    }
    return sum;
}

ll get_ans(int ord[])
{
    ll sum = 0;
    for (int i = 1; i <= m; ++i)
        for (int j = i + 1; j <= m; ++j)
            sum += c[ord[i]][ord[j]];
    return sum;
}

ll dp[N];
int fr[N], opt[N];

//对ord这个解进行下降搜索，基于Peng的dp方案
ll work(int ord[],int bg=0)
{
    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= m; ++j)
            p[i][j] = p[i][j - 1] + c[ord[i]][ord[j]];
    }
    dp[0] = 0;
    for (int i=1;i<bg;++i)
    {
        dp[i]=dp[i-1];
        fr[i]=1;
        opt[i]=0;
    }
    for (int i = bg; i <= m; ++i)
    {
        dp[i] = dp[i - 1];
        fr[i] = 1;
        opt[i] = 0;
        // cerr<<"."<<i;
        // if (i<=bg)
        //     continue;
        for (int j = 1; j <= i - 1; ++j)
        {
            //将第j个点插入到第i个点的后一个位置
            if (dp[j - 1] + (p[j][i] - p[j][j]) > dp[i])
            {
                dp[i] = dp[j - 1] + (p[j][i] - p[j][j]);
                fr[i] = i - j + 1;
                opt[i] = 1;
            }
            //将第i个点插入到第j个点的前一个位置
            if (dp[j - 1] - (p[i][i] - p[i][j - 1]) > dp[i])
            {
                dp[i] = dp[j - 1] - (p[i][i] - p[i][j - 1]);
                fr[i] = i - j + 1;
                opt[i] = 2;
            }
            ////将第i个点和第j个点交换
            if (dp[j - 1] + (p[j][i] - p[j][j]) - (p[i][i] - p[i][j]) > dp[i])
            {
                dp[i] = dp[j - 1] + (p[j][i] - p[j][j]) - (p[i][i] - p[i][j]);
                fr[i] = i - j + 1;
                opt[i] = 3;
            }

            if (j+3<i)
            {
                 //将j和j+1插入到i后
                if (dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1]);
                    fr[i]=i-j+1;
                    opt[i]=5;
                }
                 //将i-1和i插入到j前
                if (dp[j-1]-(p[i][i-2]-p[i][j-1])-(p[i-1][i-2]-p[i-1][j-1])>dp[i])
                {
                    dp[i]=dp[j-1]-(p[i][i-2]-p[i][j-1])-(p[i-1][i-2]-p[i-1][j-1]);
                    fr[i]=i-j+1;
                    opt[i]=6;
                }
            }
            if (j+4<i)
            {
                //将j,j+1和i-1,i交换
                if (dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i-2]-p[i][j+1])-(p[i-1][i-2]-p[i-1][j+1])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i-2]-p[i][j+1])-(p[i-1][i-2]-p[i-1][j+1]);
                    fr[i]=i-j+1;
                    opt[i]=7;
                }
                //将j,j+1和i交换
                if (dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i]-p[i][j+1])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i]-p[i][j+1]);
                    fr[i]=i-j+1;
                    opt[i]=8;
                }
                //将j和i-1,i交换
                if (dp[j-1]+(p[j][i]-p[j][j])-(p[i][i-2]-p[i][j])-(p[i-1][i-1]-p[i-1][j])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j])-(p[i][i-2]-p[i][j])-(p[i-1][i-1]-p[i-1][j]);
                    fr[i]=i-j+1;
                    opt[i]=9;
                }
            }
            
        }
    }
    
    if (dp[m] <= 0)
        return 0;
    int cnt = m;
    for (int i = m; i; i -= fr[i])
    {
        // cerr<<"?"<<i;
        if (fr[i] == 1)
        {
            _ord[cnt--] = ord[i];
            continue;
        }
        else
        {
            // cerr<<'('<<i-fr[i]+1<<' '<<i<<' '<<opt[i]<<' '<<dp[i-fr[i]]<<' '<<dp[i]<<')'<<endl;
            if (opt[i]<=4)
            {
                if (opt[i] & 1)
                    _ord[cnt--] = ord[i - fr[i] + 1];
                if ((opt[i] & 2) == 0)
                    _ord[cnt--] = ord[i];
                for (int j = i - 1; j > i - fr[i] + 1; --j)
                    _ord[cnt--] = ord[j];
                if ((opt[i] & 1) == 0)
                    _ord[cnt--] = ord[i - fr[i] + 1];
                if (opt[i] & 2)
                    _ord[cnt--] = ord[i];
            }
            if (opt[i]==5)
            {
                _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i;j>i-fr[i]+2;--j)
                    _ord[cnt--]=ord[j];
            }
            if (opt[i]==6)
            {
                for (int j=i-2;j>i-fr[i];--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                _ord[cnt--] = ord[i-1];
            }
            if (opt[i]==7)
            {
                _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i-2;j>i-fr[i]+2;--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                _ord[cnt--] = ord[i-1];
            }
            if (opt[i]==8)
            {
                _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i-1;j>i-fr[i]+2;--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                // _ord[cnt--] = ord[i-1];
            }
            if (opt[i]==9)
            {
                // _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i-2;j>i-fr[i]+1;--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                _ord[cnt--] = ord[i-1];
            }
        }
    }
    // cerr<<"!"<<cnt<<endl;
    for (int i = 1; i <= m; ++i)
        ord[i] = _ord[i];
    // cerr<<"!"<<dp[m]<<endl;
    return dp[m];
}
ll workbig(int ord[])
{
    dp[0] = 0;
    for (int i = 1; i <= m; ++i)
        dp[i] = -1;
    for (int i = 1; i <= m&&!tle; ++i)
    {
        // if (tle)
        //     return 0;
        // cout<<i<<endl;
        if (dp[i] < dp[i - 1])
        {
            opt[i] = 0;
            fr[i] = 1;
            dp[i] = dp[i - 1];
        }
        ll sum = 0;
        for (int j = i - 1; j&&!tle; --j)
        {
            sum -= calc(ord[i], ord[j]) - calc(ord[j], ord[i]);
            if (dp[j - 1] + sum > dp[i])
            {
                dp[i] = dp[j - 1] + sum;
                fr[i] = i - j + 1;
                opt[i] = 2;
            }
        }
        sum = 0;
        for (int j = i + 1; j <= m&&!tle; ++j)
        {
            sum += calc(ord[i], ord[j]) - calc(ord[j], ord[i]);
            if (dp[i - 1] + sum > dp[j])
            {
                dp[j] = dp[i - 1] + sum;
                fr[j] = j - i + 1;
                opt[j] = 1;
            }
        }
    }
    if (tle)
        return 0;
    if (dp[m] <= 0)
        return 0;
    int cnt = m;
    for (int i = m; i; i -= fr[i])
    {
        if (fr[i] == 1)
        {
            _ord[cnt--] = ord[i];
            continue;
        }
        else
        {
            if (opt[i] & 1)
                _ord[cnt--] = ord[i - fr[i] + 1];
            if ((opt[i] & 2) == 0)
                _ord[cnt--] = ord[i];
            for (int j = i - 1; j > i - fr[i] + 1; --j)
                _ord[cnt--] = ord[j];
            if ((opt[i] & 1) == 0)
                _ord[cnt--] = ord[i - fr[i] + 1];
            if (opt[i] & 2)
                _ord[cnt--] = ord[i];
        }
    }
    for (int i = 1; i <= m; ++i)
        ord[i] = _ord[i];
    // cout<<dp[m]<<':';for (int i=1;i<=m;++i)cout<<ord[i]<<' ';cout<<endl;
    return dp[m];
}
mt19937 rnd(0);
#ifdef exact
const int nn=20;
#else
const int nn = 18;
#endif
ll f[1 << nn];
int _fr[1 << nn];
const ll inf = 1e15;
int changeord(int ord[])
{

    int len = min(nn, m);
    int x;
    if (len == m)
        x = 0;
    else
        x = rnd() % (m - len);
    x += 1;
    ll sum = 0;
    for (int i = 0; i < len; ++i)
        for (int j = i + 1; j < len; ++j)
            sum += c[ord[i + x]][ord[j + x]];
    f[0] = 0;
    for (int i = 1; i < (1 << len); ++i)
    {
        f[i] = inf;
        for (int j = 0; j < len; ++j)
            if ((i >> j) & 1)
            {
                ll det = 0;
                for (int k = 0; k < len; ++k)
                    if ((i >> k) & 1)
                        det += c[ord[k + x]][ord[j + x]];
                if (det + f[i ^ (1 << j)] < f[i])
                {
                    _fr[i] = j;
                    f[i] = f[i ^ (1 << j)] + det;
                }
            }
    }
    if (f[(1 << len) - 1] == sum)
        return 0;
    vector<int> nord;
    for (int i = (1 << len) - 1; i;)
    {
        nord.push_back(ord[_fr[i] + x]);
        i -= (1 << _fr[i]);
    }
    assert(nord.size() == len);
    for (int i = 0; i < len; ++i)
        ord[i + x] = nord[len - i - 1];

    // set<int> ss;
    // for (int i=0;i<n)
    return x;
}
int cc[nn][nn];
int changeordbig(int ord[])
{

    int len = min(nn, m);
    int x;
    if (len == m)
        x = 0;
    else
        x = rnd() % (m - len);
    x += 1;
    ll sum = 0;
    for (int i=0;i<len;++i)
        for (int j=0;j<len;++j)
            if (i!=j)
                cc[i][j]=calc(ord[x+i],ord[x+j]);
    for (int i = 0; i < len; ++i)
        for (int j = i + 1; j < len; ++j)
            sum +=cc[i][j];// c[ord[i + x]][ord[j + x]];
    f[0] = 0;
    for (int i = 1; i < (1 << len); ++i)
    {
        f[i] = inf;
        for (int j = 0; j < len; ++j)
            if ((i >> j) & 1)
            {
                ll det = 0;
                for (int k = 0; k < len; ++k)
                    if ((i >> k) & 1)
                        det += cc[k][j];//c[ord[k + x]][ord[j + x]];
                if (det + f[i ^ (1 << j)] < f[i])
                {
                    _fr[i] = j;
                    f[i] = f[i ^ (1 << j)] + det;
                }
            }
    }
    if (f[(1 << len) - 1] == sum)
        return 0;
    vector<int> nord;
    for (int i = (1 << len) - 1; i;)
    {
        nord.push_back(ord[_fr[i] + x]);
        i -= (1 << _fr[i]);
    }
    assert(nord.size() == len);
    for (int i = 0; i < len; ++i)
        ord[i + x] = nord[len - i - 1];

    // set<int> ss;
    // for (int i=0;i<n)
    return x;
}
int swapsum=0;
void cross(int parent1[], int parent2[], int child[])
{
    int x = rnd() % m + 1;
    int y = rnd() % m + 1;
    if (x > y)
        swap(x, y);
    for (int i = 0; i < m; ++i)
        used[i] = 0;
    for (int i = x; i <= y; ++i)
    {
        child[i] = parent1[i];
        used[child[i]] = 1;
    }
    int pos = 0;
    for (int i = 1; i <= m; ++i)
    {
        if (!used[parent2[i]])
        {
            ++pos;
            if (pos == x)
                pos = y + 1;
            child[pos] = parent2[i];
        }
    }
    // for (int i=1;i<=swapsum;++i)
    // {
    //     int x=rnd()%m+1;int y=rnd()%m+1;
    //     if (x!=y)
    //         swap(child[x],child[y]);
    // }
    // swapsum+=5;
    // if (swapsum>n/5)
    //     swapsum=0;
}
ll get_hash(int ord[])
{
    ll sum = 0;
    for (int i = 1; i <= m; ++i)
        sum = ((ll)sum * base % mod + ord[i]) % mod;
    return sum;
}

int main(int argc, char *argv[])
{
    // freopen("1.gr","r",stdin);
#ifdef sig
    struct sigaction action; // 信号
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);
#endif
    clock_t clockBegin = clock();
    read(n);                          // A点数
    read(m);                          // B点数
    read(K);                          // 边数
    for (int i = 0, u, v; i < K; ++i) // 每个点从0开始编号
    {
        read(u);
        --u;
        read(v);
        --v;
        v -= n;
        e[v].push_back(u); // e保存了所有的边，e[v]表示v的邻居 v in B
    }
    for (int i = 0; i < m; ++i) // 遍历b
        if (!e[i].empty())
        {
            sort(e[i].begin(), e[i].end()); // i的邻居升序排序
            for (int x : e[i])
                avg[i] += x;
            avg[i] /= (double)e[i].size(); // avg[i]是i的邻居的平均数
        }
        else
        {
            avg[i] = 0;
            // cerr<<"!"<<i<<endl;
        }
    for (int i = 1; i <= m; ++i)
        ord[i] = i - 1;
    sort(ord + 1, ord + 1 + m, cmp); // B的顺序，初始解； cmp比较的是avg
    // cerr<<m<<' '<<M<<endl;
    /**
    将B>21000时，认为内存存不下；N^2空间
    */
    if (m < M)
    // if (0)
    {
        for (int i = 0; i < m; ++i)
        {
            c[i].resize(m);
            p[i].resize(m + 1);
        }
        p[m].resize(m + 1);
        for (int i = 0; i < m; ++i)
            for (int j = i + 1; j < m; ++j)
            {
                ll pi = calc(i, j), pj = calc(j, i);
                c[i][j] = pi - pj; // C[i][j] 点i在j之前的，N[i]和N[j]产生的交叉数
                c[j][i] = pj - pi; //
            }
        // cerr<<"finish"<<' '<<clock()-clockBegin<<endl;
#ifdef exact
int mxt=29*60;
#else
int mxt = 280;   
#endif  
          // 29*60; //时间cutoff
        // cerr<<60*30*CLOCKS_PER_SEC<<endl;
        int mx = 0;                  // 种群大小，不是常数，会变化；最大10
        for (int i=0;i==0||(i<10&&clock()-clockBegin<=mxt/5*CLOCKS_PER_SEC);++i) // 最大的种群为10，构造10次初始解
        {
            //  cerr << i << endl;
            for (int j = 1; j <= m; ++j)
                ords[i][j] = ord[j]; // orders[i]表示第i个个体
            if (i >= 1) //除了orders[0]，其他的解都随机化
            {
                shuffle(ords[i] + 1, ords[i] + 1 + m, rnd);
                // for (int j=1;j<=m;++j)
                //     swap(ords[i][j],ords[i][j+rnd()%min(m-j+1,i*10)]);
            }

            ll x = get_hash(ords[i]); 
            int cnt = 0;
            /**
             * tabu.count(x) != 0 ->x出现过            
            */
            while (tabu.count(x) != 0 && cnt <= 10 && !tle && clock() - clockBegin <= mxt / 5 * CLOCKS_PER_SEC)
            {
                shuffle(ords[i] + 1, ords[i] + 1 + m, rnd);
                x = get_hash(ords[i]);
                ++cnt;
            }
            // cerr<<"hashb"<<x<<endl;
            //结束构造初始种群
            if (tabu.count(x) != 0 || tle)// || clock() - clockBegin > mxt / 2 * CLOCKS_PER_SEC)
                break;
            // tabu.insert(x);
            // ll tmp;
            // while ((tmp=work(ords[i]))&& !tle)// && clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC)
            //     cerr<<tmp<<' '<<get_ans(ords[i])<<' '<<(clock()-clockBegin)/(double)CLOCKS_PER_SEC<<';';

            while (work(ords[i]) && !tle && clock() - clockBegin <= mxt / 5 * CLOCKS_PER_SEC)
                ;
            x = get_hash(ords[i]);
            while (tabu.count(x) != 0 && cnt <= 10 && !tle && clock() - clockBegin <= mxt / 5 * CLOCKS_PER_SEC)
            {
                shuffle(ords[i] + 1, ords[i] + 1 + m, rnd);
                while (work(ords[i]) && !tle && clock() - clockBegin <= mxt / 5 * CLOCKS_PER_SEC)
                    ;
                x = get_hash(ords[i]);
                ++cnt;
            }
            if (tle)
                break;
            if (tabu.count(x))
                break;
            
            tabu.insert(x);
            if (!tle)
            {
                mx = i + 1;
                ans[i] = get_ans(ords[i]);
            }
            // cerr<<"hash:"<<x<<' '<<ans[i]<<endl;
            // cerr<<(clock()-clockBegin)/(double)CLOCKS_PER_SEC<<endl;
        }
        // cerr<<"!"<<clock()-clockBegin<<endl;
        if (mx == 0)
        {
            for (int i = 1; i <= m; ++i)
                printf("%d\n", ord[i] + n + 1);
            return 0;
        }
        int pc=0;
        while (clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC && !tle && mx > 1)
        {
            // cerr<<(clock()-clockBegin)/(double)CLOCKS_PER_SEC<<endl;
            // cerr<<'.';
            int x = rnd() % mx, y = rnd() % mx;
            while (x == y)
                y = rnd() % mx;
            cross(ords[x], ords[y], child);
            int z = get_hash(child);
            if (tabu.count(z))
                continue;
            while (work(child) && !tle && clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC)
                ;
            ll tmphash = get_hash(child);
            if (tabu.count(tmphash))
                continue;
            int mnp = 0;
            tabu.insert(tmphash);
            for (int i = 1; i < mx; ++i)
                if (ans[i] > ans[mnp])// || (ans[i] == ans[mnp] && rnd() % 2 == 0))
                    mnp = i;
            ll tmpans = get_ans(child);
            // if (++pc<=10)cerr<<tmphash<<' '<<tmpans<<endl;
            if (tmpans <= ans[mnp])
            {
                ans[mnp] = tmpans;
                for (int i = 1; i <= m; ++i)
                    ords[mnp][i] = child[i];
                // cout<<tmpans<<' '<<tmphash<<' ';for (int i=1;i<=m;++i)cout<<child[i]<<' ';cout<<endl;
            }
        }
        int mxp = 0;
        for (int i = 1; i < mx; ++i)
            if (ans[i] < ans[mxp])
                mxp = i;
        for (int i = 1; i <= m; ++i)
            ord[i] = ords[mxp][i];
        // cerr<<"!"<<get_ans(ord)<<endl;
        // while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
        // {

        //     if (!work(ord))
        //     {
        //         while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
        //         {
        //             if (changeord(ord))
        //                 break;
        //         }
        //     }
        //     // break;
        // }
        int bg=0;
        while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
        {
            // cerr<<clock()-clockBegin<<endl;
            if (!work(ord,bg))
            {
                while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
                {
                    int pos=changeord(ord);
                    if (pos!=0)
                    {
                        bg=pos;
                        break;
                    }
                }
            }
            else
                bg=0;
            // break;
        }
    }
    else
    {
#ifdef exact
        int mxt=20*60;
        while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
        {
            if (!workbig(ord))
            {
                while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
                {
                    if (changeordbig(ord))
                        break;
                }
            }
            // break;
        }
#else
        while (clock() - clockBegin <= 5 * CLOCKS_PER_SEC && !tle)
        {
            if (!workbig(ord))
                break;
        }
#endif
    }
    for (int i = 1; i <= m; ++i)
        printf("%d\n", ord[i] + n + 1);
    // printf("%lld\n",get_ans(ord));
    cerr<<"!"<<get_ans(ord)<<endl;
    return 0;
}