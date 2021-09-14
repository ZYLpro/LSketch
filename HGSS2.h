#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <cmath>
#include <stdlib.h>
#include <bitset>
#include <memory.h>
#include <time.h>
#include <ctime>
#ifndef HASHTABLE_H
#define HASHTABLE_H
#include "head/GraphUndirected.h"
#endif // HASHTABLE_H
using namespace std;
#define prime 739
#define bigger_p 1048576
#define timer 5
#define M 80000

#define Roomnum 2   // This is the parameter to controll the maximum number of rooms in a bucket.
#define sub_num 168 // subwindow num

struct cell
{
    unsigned short src[Roomnum];
    unsigned short dst[Roomnum];
    unsigned int idx;
    short C[Roomnum][sub_num]; // stores the total count
    int P[Roomnum][sub_num];   // stores the edge prime
    unsigned short h1[Roomnum];
    unsigned short h2[Roomnum];
};
struct mnode
{
    unsigned int h;
    unsigned short g;
};
struct link_node
{
    unsigned int key;
    short C[sub_num] = {};
    int P[sub_num] = {};
    int label = 1;      //
    time_t last_T = -1; //stores the starting time for last sub window
    link_node *next;
};

class HGSS
{
private:
    hashTable<string> mapTable;
    int w;
    int b; // block
    int r;
    int p;
    int s;
    int f;
    int sub_size; // sub window size
    int win_k;    // sub window num
    bool useT;
    int tablesize;

    cell *value;
    time_t lastT;             //stores the starting time for last sub window
    map<int, string> k2label; // not used

public:
    int prime_num = 168; //control #of prime numbers we use
    int Prime[168] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                      31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                      73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
                      127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
                      179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
                      233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
                      283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
                      353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
                      419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
                      467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
                      547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
                      607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
                      661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
                      739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
                      811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
                      877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
                      947, 953, 967, 971, 977, 983, 991, 997};
    map<string, int> random_index; 

    vector<link_node *> buffer;
    map<unsigned int, int> index;
    int buffer_index;
    int edge_num; // count the number of edges in the buffer to assist buffer size analysis. Self loop edge is not included as it does not use additional memory.
    HGSS(int width, int block, int range, int p_num, int size, int f_num, int sub_s, int window_k, bool usetable, int tablesize = 0);
    ~HGSS()
    {
        delete[] value;
        cleanupBuffer();
    }
    // insert
    void insert(string s1, string s2, int weight, string l_e, time_t time);                         // edge_label   road
    void insert(string s1, string s2, int weight, string l_A, string l_B, string l_e, time_t time); // node_label1 node_label2 edge_label
    void cleanupBuffer();

    // edge query
    // road
    int edgeQuery(string s1, string s2, time_t end_time);
    int edgeQueryWithLabel(string s1, string s2, string label, time_t end_time);
    // phone
    int edgeQuery(string s1, string l_A, string s2, string l_B, time_t end_time);
    int edgeQueryWithLabel(string s1, string l_A, string s2, string l_B, string label, time_t end_time);

    // pathReachability
    // road
    bool query(string s1, string s2, time_t end_time);
    bool queryWithLabel(string s1, string s2, string l_e, time_t end_time);
    // phone  unfinished
    bool query(string s1, string l_A, string s2, string l_B, time_t end_time);
    bool queryWithLabel(string s1, string l_A, string s2, string l_B, string l_e, time_t end_time);

    //  node out
    int nodeValueQuery(string s1, int type, time_t end_time);             //src_type = 0 dst_type = 1
    int nodeValueQuery(string s1, int type, string l_A, time_t end_time); //src_type = 0 dst_type = 1
    // 
    int nodeValueQueryWithLabel(string s1, string l_e, int type, time_t end_time);             //src_type = 0 dst_type = 1
    int nodeValueQueryWithLabel(string s1, string l_e, int type, string l_A, time_t end_time); //src_type = 0 dst_type = 1

    int nodeDegreeQuery(string s1, int type); //src_type = 0 dst_type = 1
    void nodeSuccessorQuery(string s1, vector<string> &IDs);
    void nodePrecursorQuery(string s2, vector<string> &IDs);
    int TriangleCounting();
};

HGSS::HGSS(int width, int block, int range, int p_num, int size, int f_num, int sub_s, int window_k, bool usehashtable, int TableSize) //the side length of matrix, the length of hash addtress list, the number of candidate bucekt
// the number of rooms, whether to use hash table, and the size of the table.
// Hash table which stores the original nodes can be omitted if not needed. For nodequery,
//  reachability, edgequery not needed. But needed for triangel counting, degree query, and successor / precursor queries.
{
    w = width;
    b = block;
    r = range;        /* r x r mapped baskets */
    p = p_num;        /*candidate buckets*/
    s = size;         /*multiple rooms*/
    f = f_num;        /*finger print lenth*/
    sub_size = sub_s; //  subwindow size
    win_k = window_k; //  subwindow num
    buffer_index = 0;
    edge_num = 0;
    value = new cell[w * w];
    lastT = -1;
    useT = usehashtable;
    tablesize = TableSize;
    memset(value, 0, sizeof(cell) * w * w);
    if (usehashtable)
        mapTable.init(tablesize);
}

void HGSS::cleanupBuffer()
{
    vector<link_node *>::iterator IT = buffer.begin();
    link_node *e, *tmp;
    for (; IT != buffer.end(); ++IT)
    {
        e = *IT;
        while (e != NULL)
        {
            tmp = e->next;
            delete e;
            e = tmp;
        }
    }
}
void HGSS::insert(string s1, string s2, int weight, string l_e, time_t time) // road, without node labels
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    unsigned int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;

    unsigned int k1 = (h1 << f) + g1;
    unsigned int k2 = (h2 << f) + g2;

    if (useT)
    {
        mapTable.insert(k1, s1);
        mapTable.insert(k2, s2);
    }

   
    int prime_group_num = 5; 
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    map<string, int>::iterator it0 = random_index.find(l_e);
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    int *tmp1 = new int[r];
    int *tmp2 = new int[r];
    tmp1[0] = g1;
    tmp2[0] = g2;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
        tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
    }
    bool inserted = false;
    long key = g1 + g2;

    // sliding windows
    if (lastT > 0)
    {
        double cost = difftime(time, lastT);
        while ((int)cost >= sub_size) //value[pos].lastT + sub_s ==  time
        {
            for (int i = 0; i < w * w; i++)
            {
                for (int j = 0; j < s; j++)
                {
                    if (value[i].src[j] != 0 && value[i].dst[j] != 0)  //slides
                    {
                        for (int u = 1; u < win_k; u++)
                        {
                            value[i].C[j][u - 1] = value[i].C[j][u];
                            value[i].P[j][u - 1] = value[i].P[j][u];
                        }
                        value[i].C[j][win_k - 1] = 0;
                        value[i].P[j][win_k - 1] = 1; 
                    }
                }
            }
            lastT = lastT + sub_size;
            cost -= sub_size;
        }
    }
    else
        lastT = time; 

    for (int i = 0; i < p; i++)
    {
        key = (key * timer + prime) % bigger_p;
        int index = key % (r * r);
        int index1 = index / r;
        int index2 = index % r;
        int p1 = (h1 + tmp1[index1]) % w;
        int p2 = (h2 + tmp2[index2]) % w;

        int pos = p1 * w + p2;
        for (int j = 0; j < s; j++)
        {
            if ((((value[pos].idx >> (j << 3)) & ((1 << 8) - 1)) == (index1 | (index2 << 4))) && (value[pos].src[j] == g1) && (value[pos].dst[j] == g2))
            {
                value[pos].C[j][win_k - 1] += weight;
                value[pos].P[j][win_k - 1] *= edge_label;
                inserted = true;
                break;
            }
            if (value[pos].src[j] == 0)
            {
                value[pos].idx |= ((index1 | (index2 << 4)) << (j << 3));
                value[pos].src[j] = g1;
                value[pos].dst[j] = g2;
                // h
                value[pos].h1[j] = h1;
                value[pos].h2[j] = h2;

                value[pos].C[j][win_k - 1] = weight;
                value[pos].P[j][win_k - 1] = edge_label;

                inserted = true;
                break;
            }
        }
        if (inserted)
            break;
    }
    if (!inserted)
    {
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];  
            node->label *= edge_label; 
            while (true)
            {
                if (node->key == k2)
                {
                    if (node->last_T != -1) 
                    {
                        double cost = difftime(time, node->last_T);
                        while ((int)cost >= sub_size)
                        {
                            if (node->P[0] != 0)
                                node->label = node->label / node->P[0]; 
                            for (int u = 1; u < win_k; u++)
                            {
                                node->C[u - 1] = node->C[u];
                                node->P[u - 1] = node->P[u];
                            }
                            node->C[win_k - 1] = 0;
                            node->P[win_k - 1] = 1;

                            node->last_T = node->last_T + sub_size;
                            cost -= sub_size;
                        }
                        node->C[win_k - 1] += weight;
                        node->P[win_k - 1] *= edge_label;
                    }
                    else
                    {
                        node->last_T = time;
                        node->label = edge_label;

                        node->C[win_k - 1] = weight;
                        node->P[win_k - 1] = edge_label;
                    }
                    break;
                }
                if (node->next == NULL)
                {
                    link_node *ins = new link_node;
                    ins->key = k2;
                    ins->C[win_k - 1] = weight;
                    ins->P[win_k - 1] = edge_label;
                    ins->last_T = time;
                    ins->next = NULL;
                    node->next = ins;
                    edge_num++;
                    break;
                }
                node = node->next;
            }
        }
        else
        {
            index[k1] = buffer_index;
            buffer_index++;
            link_node *node = new link_node;
            node->key = k1;
            node->label = edge_label;
            if (k1 != k2) //k1==k2 means loop
            {
                link_node *ins = new link_node;
                ins->key = k2;
                ins->C[win_k - 1] = weight;
                ins->P[win_k - 1] = edge_label;
                ins->last_T = time;
                ins->next = NULL;
                node->next = ins;
                edge_num++;
            }
            else
            {
                node->C[win_k - 1] = weight;
                node->P[win_k - 1] = edge_label;
                node->last_T = time; 
                node->next = NULL;
            }
            buffer.push_back(node);
        }
    }
    delete[] tmp1;
    delete[] tmp2;
    return;
}

void HGSS::insert(string s1, string s2, int weight, string l_A, string l_B, string l_e, time_t time)
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    unsigned int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;

    unsigned int k1 = (h1 << f) + g1;
    unsigned int k2 = (h2 << f) + g2;

    if (useT)
    {
        mapTable.insert(k1, s1);
        mapTable.insert(k2, s2);
    }

    // blocks
    unsigned int label1 = (*hfunc[0])((unsigned char *)(l_A.c_str()), l_A.length());
    unsigned int label2 = (*hfunc[0])((unsigned char *)(l_B.c_str()), l_B.length());
    unsigned int m = (label1 % (w / b));
    unsigned int n = (label2 % (w / b));

    
    auto iter1 = k2label.find(k1);
    if (iter1 == k2label.end())
        k2label[k1] = l_A;
    auto iter2 = k2label.find(k2);
    if (iter2 == k2label.end())
        k2label[k2] = l_B;
    int prime_group_num = 5; 
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    map<string, int>::iterator it0 = random_index.find(l_e);
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    int *tmp1 = new int[r];
    int *tmp2 = new int[r];
    tmp1[0] = g1;
    tmp2[0] = g2;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
        tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
    }
    bool inserted = false;
    long key = g1 + g2;

    // sliding windows
    if (lastT > 0)
    {
        double cost = difftime(time, lastT);
        while ((int)cost >= sub_size) //value[pos].lastT + sub_s ==  time
        {
            for (int i = 0; i < w * w; i++)
            {
                for (int j = 0; j < s; j++)
                {
                    if (value[i].src[j] != 0 && value[i].dst[j] != 0)
                    {
                        for (int u = 1; u < win_k; u++)
                        {
                            value[i].C[j][u - 1] = value[i].C[j][u];
                            value[i].P[j][u - 1] = value[i].P[j][u];
                        }
                        value[i].C[j][win_k - 1] = 0;
                        value[i].P[j][win_k - 1] = 1; // 
                    }
                }
            }
            lastT = lastT + sub_size;
            cost -= sub_size;
        }
    }
    else
        lastT = time; 

    for (int i = 0; i < p; i++)
    {
        key = (key * timer + prime) % bigger_p;
        int index = key % (r * r);
        int index1 = index / r;
        int index2 = index % r;
        int p1 = (h1 + tmp1[index1]) % b; //在block中的位置
        int p2 = (h2 + tmp2[index2]) % b;

        int pos = (m * b + p1) * w + (n * b + p2);
        for (int j = 0; j < s; j++)
        {
            if ((((value[pos].idx >> (j << 3)) & ((1 << 8) - 1)) == (index1 | (index2 << 4))) && (value[pos].src[j] == g1) && (value[pos].dst[j] == g2))
            {
                value[pos].C[j][win_k - 1] += weight;
                value[pos].P[j][win_k - 1] *= edge_label;
                inserted = true;
                break;
            }
            if (value[pos].src[j] == 0)
            {
                value[pos].idx |= ((index1 | (index2 << 4)) << (j << 3));
                value[pos].src[j] = g1;
                value[pos].dst[j] = g2;
                // h
                value[pos].h1[j] = h1;
                value[pos].h2[j] = h2;

                value[pos].C[j][win_k - 1] = weight;
                value[pos].P[j][win_k - 1] = edge_label;

                inserted = true;
                break;
            }
        }
        if (inserted)
            break;
    }
    if (!inserted)
    {
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            link_node *head = buffer[tag];
            head->label *= edge_label; 
            while (true)
            {
                if (node->key == k2)
                {
                    if (node->last_T != -1) 
                    {
                        double cost = difftime(time, node->last_T);
                        while ((int)cost >= sub_size)
                        {
                            if (node->P[0] != 0)
                                head->label = head->label / node->P[0]; 
                            for (int u = 1; u < win_k; u++)
                            {
                                node->C[u - 1] = node->C[u];
                                node->P[u - 1] = node->P[u];
                            }
                            node->C[win_k - 1] = 0;
                            node->P[win_k - 1] = 1;

                            node->last_T = node->last_T + sub_size;
                            cost -= sub_size;
                        }
                        node->C[win_k - 1] += weight;
                        node->P[win_k - 1] *= edge_label;
                    }
                    else
                    {
                        node->last_T = time;
                        node->label *= edge_label;

                        node->C[win_k - 1] = weight;
                        node->P[win_k - 1] = edge_label;
                    }
                    break;
                }
                if (node->next == NULL)
                {
                    link_node *ins = new link_node;
                    ins->key = k2;
                    ins->C[win_k - 1] = weight;
                    ins->P[win_k - 1] = edge_label;
                    ins->last_T = time;
                    ins->next = NULL;
                    node->next = ins;
                    edge_num++;
                    break;
                }
                node = node->next;
            }
        }
        else
        {
            index[k1] = buffer_index;
            buffer_index++;
            link_node *node = new link_node;
            node->key = k1;
            node->label = edge_label;
            if (k1 != k2) //k1==k2 means loop
            {
                link_node *ins = new link_node;
                ins->key = k2;
                ins->C[win_k - 1] = weight;
                ins->P[win_k - 1] = edge_label;
                ins->last_T = time;
                ins->next = NULL;
                node->next = ins;
                edge_num++;
            }
            else 
            {
                node->C[win_k - 1] = weight;
                node->P[win_k - 1] = edge_label;
                node->last_T = time; 
                node->next = NULL;
            }
            buffer.push_back(node);
        }
    }
    delete[] tmp1;
    delete[] tmp2;
    return;
}
/*type 0 is for successor query, type 1 is for precusor query*/
int HGSS::nodeValueQuery(string s1, int type, time_t end_time) // road
{
    int weight = 0;
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    int *tmp1 = new int[r];
    tmp1[0] = g1;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
    }
    for (int i = 0; i < r; i++)
    {
        int p1 = (h1 + tmp1[i]) % w;
        for (int k = 0; k < w; k++)
        {
            if (type == 0) /*successor query*/
            {
                int pos = p1 * w + k;
                for (int j = 0; j < s; ++j)
                {
                    if (type == 0 && (((value[pos].idx >> ((j << 3))) & ((1 << 4) - 1)) == i) && (value[pos].src[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                            weight += value[pos].C[j][t];
                    }
                }
            }
            else if (type == 1) /*precursor query*/
            {
                int pos = p1 + k * w;
                for (int j = 0; j < s; ++j)
                {
                    if (type == 1 && (((value[pos].idx >> ((j << 3) + 4)) & ((1 << 4) - 1)) == i) && (value[pos].dst[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                            weight += value[pos].C[j][t];
                    }
                }
            }
        }
    }
    if (type == 0)
    {
        unsigned int k1 = (h1 << f) + g1;
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
                weight += node->C[t];
            node = node->next;
            while (node != NULL)
            {
                int b_i2 = floor(difftime(end_time, node->last_T) / sub_size);
                for (int t = b_i2; t < win_k; t++)
                    weight += node->C[t];
                node = node->next;
            }
        }
    }
    else if (type == 1)
    {
        unsigned int k1 = (h1 << f) + g1;
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int b_i1 = floor(difftime(end_time, buffer[it->second]->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
                weight += buffer[it->second]->C[t];
        }
        // others
        for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            node = node->next;
            while (node != NULL)
            {
                if (node->key == k1)
                {
                    int b_i2 = floor(difftime(end_time, node->last_T) / sub_size);
                    for (int t = b_i2; t < win_k; t++)
                        weight += node->C[t];
                }
                node = node->next;
            }
        }
    }
    delete[] tmp1;
    return weight;
}

int HGSS::nodeValueQueryWithLabel(string s1, string l_e, int type, time_t end_time) // road
{
    int weight = 0;
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    int *tmp1 = new int[r];
    tmp1[0] = g1;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
    }

    // 引入质数分组
    int prime_group_num = 5; //一个组包含多少个质数
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    // 第一次出现，随机
    // 第二次及以后，使用第一次的随机数（否则同一个边label会对应不同的质数）
    map<string, int>::iterator it0 = random_index.find(l_e);
    // 存在label的随机数
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    // 最后一组可能会越界
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    for (int i = 0; i < r; i++)
    {
        int p1 = (h1 + tmp1[i]) % w;
        for (int k = 0; k < w; k++)
        {
            if (type == 0) /*successor query*/
            {
                int pos = p1 * w + k;
                for (int j = 0; j < s; ++j)
                {
                    if (type == 0 && (((value[pos].idx >> ((j << 3))) & ((1 << 4) - 1)) == i) && (value[pos].src[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[j][t];
                            while (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                weight++;
                                temp_P = temp_P / edge_label;
                            }
                        }
                    }
                }
            }
            else if (type == 1) /*precursor query*/
            {
                int pos = p1 + k * w;
                for (int j = 0; j < s; ++j)
                {
                    if (type == 1 && (((value[pos].idx >> ((j << 3) + 4)) & ((1 << 4) - 1)) == i) && (value[pos].dst[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[j][t];
                            while (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                weight++;
                                temp_P = temp_P / edge_label;
                            }
                        }
                    }
                }
            }
        }
    }
    if (type == 0)
    {
        unsigned int k1 = (h1 << f) + g1;
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            // 自连边的重量
            int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
            {
                int temp_P = node->P[t];
                while (temp_P != 0 && temp_P % edge_label == 0)
                {
                    weight++;
                    temp_P = temp_P / edge_label;
                }
            }
            node = node->next;
            while (node != NULL)
            {
                int b_i2 = floor(difftime(end_time, node->last_T) / sub_size);
                for (int t = b_i2; t < win_k; t++)
                {
                    int temp_P = node->P[t];
                    while (temp_P != 0 && temp_P % edge_label == 0)
                    {
                        weight++;
                        temp_P = temp_P / edge_label;
                    }
                }
                node = node->next;
            }
        }
    }
    else if (type == 1)
    {
        unsigned int k1 = (h1 << f) + g1;
        // 计算自连边的入度贡献
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int b_i1 = floor(difftime(end_time, buffer[it->second]->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
            {
                int temp_P = buffer[it->second]->P[t];
                while (temp_P != 0 && temp_P % edge_label == 0)
                {
                    weight++;
                    temp_P = temp_P / edge_label;
                }
            }
        }
        // 其他
        for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            node = node->next;
            while (node != NULL)
            {
                if (node->key == k1)
                {
                    int b_i2 = floor(difftime(end_time, node->last_T) / sub_size);
                    for (int t = b_i2; t < win_k; t++)
                    {
                        int temp_P = node->P[t];
                        while (temp_P != 0 && temp_P % edge_label == 0)
                        {
                            weight++;
                            temp_P = temp_P / edge_label;
                        }
                    }
                }
                node = node->next;
            }
        }
    }
    delete[] tmp1;
    return weight;
}

/*type 0 is for successor query, type 1 is for precusor query*/
int HGSS::nodeValueQuery(string s1, int type, string l_A, time_t end_time)
{
    int weight = 0;
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    int *tmp1 = new int[r];
    tmp1[0] = g1;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
    }
    // 矩阵块
    unsigned int label1 = (*hfunc[0])((unsigned char *)(l_A.c_str()), l_A.length());
    unsigned int m = (label1 % (w / b));

    for (int i = 0; i < r; i++)
    {
        int p1 = (h1 + tmp1[i]) % b; //block
        for (int k = 0; k < w; k++)
        {
            if (type == 0) /*successor query*/
            {
                int pos = (m * b + p1) * w + k;
                for (int j = 0; j < s; ++j)
                {
                    if (type == 0 && (((value[pos].idx >> ((j << 3))) & ((1 << 4) - 1)) == i) && (value[pos].src[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                            weight += value[pos].C[j][t];
                    }
                }
            }
            else if (type == 1) /*precursor query*/
            {
                int pos = k * w + (m * b + p1);
                for (int j = 0; j < s; ++j)
                {
                    if (type == 1 && (((value[pos].idx >> ((j << 3) + 4)) & ((1 << 4) - 1)) == i) && (value[pos].dst[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                            weight += value[pos].C[j][t];
                    }
                }
            }
        }
    }
    if (type == 0)
    {
        unsigned int k1 = (h1 << f) + g1;
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            // 自连边的重量
            int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
                weight += node->C[t];
            node = node->next;
            while (node != NULL)
            {
                int b_i2 = floor(difftime(end_time, node->last_T) / sub_size);
                for (int t = b_i2; t < win_k; t++)
                    weight += node->C[t];
                node = node->next;
            }
        }
    }
    else if (type == 1)
    {
        unsigned int k1 = (h1 << f) + g1;
        // 计算自连边的入度贡献
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int b_i1 = floor(difftime(end_time, buffer[it->second]->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
                weight += buffer[it->second]->C[t];
        }
        // 其他
        for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            node = node->next;
            while (node != NULL)
            {
                if (node->key == k1)
                {
                    int b_i2 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
                    for (int t = b_i2; t < win_k; t++)
                        weight += node->C[t];
                }
                node = node->next;
            }
        }
    }
    delete[] tmp1;
    return weight;
}

int HGSS::nodeValueQueryWithLabel(string s1, string l_e, int type, string l_A, time_t end_time)
{
    int weight = 0;
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    int *tmp1 = new int[r];
    tmp1[0] = g1;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
    }
    // 矩阵块
    unsigned int label1 = (*hfunc[0])((unsigned char *)(l_A.c_str()), l_A.length());
    unsigned int m = (label1 % (w / b));

    // 引入质数分组
    int prime_group_num = 5; //一个组包含多少个质数
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    // 第一次出现，随机
    // 第二次及以后，使用第一次的随机数（否则同一个边label会对应不同的质数）
    map<string, int>::iterator it0 = random_index.find(l_e);
    // 存在label的随机数
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    // 最后一组可能会越界
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    for (int i = 0; i < r; i++)
    {
        int p1 = (h1 + tmp1[i]) % b; //block
        for (int k = 0; k < w; k++)
        {
            if (type == 0) /*successor query*/
            {
                int pos = (m * b + p1) * w + k;
                for (int j = 0; j < s; ++j)
                {
                    if (type == 0 && (((value[pos].idx >> ((j << 3))) & ((1 << 4) - 1)) == i) && (value[pos].src[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[j][t];
                            while (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                weight++;
                                temp_P = temp_P / edge_label;
                            }
                        }
                    }
                }
            }
            else if (type == 1) /*precursor query*/
            {
                int pos = k * w + (m * b + p1);
                for (int j = 0; j < s; ++j)
                {
                    if (type == 1 && (((value[pos].idx >> ((j << 3) + 4)) & ((1 << 4) - 1)) == i) && (value[pos].dst[j] == g1))
                    {
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[j][t];
                            while (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                weight++;
                                temp_P = temp_P / edge_label;
                            }
                        }
                    }
                }
            }
        }
    }
    if (type == 0)
    {
        unsigned int k1 = (h1 << f) + g1;
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            // int temp = node->label;
            // while (temp != 0 && temp % edge_label == 0)
            // {
            // 	weight++;
            // 	temp = temp / edge_label;
            // }
            // 自连边的重量
            int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
            {
                int temp_P = node->P[t];
                while (temp_P != 0 && temp_P % edge_label == 0)
                {
                    weight++;
                    temp_P = temp_P / edge_label;
                }
            }
            node = node->next;
            while (node != NULL)
            {
                int b_i2 = floor(difftime(end_time, node->last_T) / sub_size);
                for (int t = b_i2; t < win_k; t++)
                {
                    int temp_P = node->P[t];
                    while (temp_P != 0 && temp_P % edge_label == 0)
                    {
                        weight++;
                        temp_P = temp_P / edge_label;
                    }
                }
                node = node->next;
            }
        }
    }
    else if (type == 1)
    {
        unsigned int k1 = (h1 << f) + g1;
        // 计算自连边的入度贡献
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int b_i1 = floor(difftime(end_time, buffer[it->second]->last_T) / sub_size); // 满足时间的窗口的起始index
            for (int t = b_i1; t < win_k; t++)
            {
                int temp_P = buffer[it->second]->P[t];
                while (temp_P != 0 && temp_P % edge_label == 0)
                {
                    weight++;
                    temp_P = temp_P / edge_label;
                }
            }
        }
        // 其他
        for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            node = node->next;
            while (node != NULL)
            {
                if (node->key == k1)
                {
                    int b_i2 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
                    for (int t = b_i2; t < win_k; t++)
                    {
                        int temp_P = node->P[t];
                        while (temp_P != 0 && temp_P % edge_label == 0)
                        {
                            weight++;
                            temp_P = temp_P / edge_label;
                        }
                    }
                }
                node = node->next;
            }
        }
    }
    delete[] tmp1;
    return weight;
}

// void HGSS::nodeSuccessorQuery(string s1, vector<string>&IDs)// query the successors of a node, s1 is the ID of the queried node. results are put in the vector, hash table needed.
// {
// 	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
// 	int tmp=pow(2,f)-1;
// 	unsigned short g1=hash1 & tmp;
// 	if(g1==0) g1+=1;
// 	unsigned int h1 = (hash1>>f)%w;
// 	unsigned int k1 = (h1 << f) + g1;
// 	int* tmp1 = new int[r];
// 	tmp1[0] = g1;
// 	for (int i = 1; i < r; i++)
// 	{
// 		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
// 	}
// 	for (int i = 0; i < r; i++)
// 	{
// 		int p1 = (h1 + tmp1[i]) % w;
// 		for (int k = 0; k < w; k++)
// 		{
// 			int pos = p1*w + k;
// 			for (int j = 0; j < s; ++j)
// 			{
// 				if ((((value[pos].idx >> ((j << 3)))&((1 << 4) - 1)) == i) && (value[pos].src[j] == g1))
// 				{
// 					     int tmp_g = value[pos].dst[j];
// 						 int tmp_s = ((value[pos].idx >> ((j << 3) + 4))&((1 << 4) - 1));

// 						 int shifter = tmp_g;
// 						 for (int v = 0; v < tmp_s; v++)
// 							 shifter = (shifter*timer + prime) % bigger_p;
// 						 int tmp_h = k;
// 						 while (tmp_h < shifter)
// 							 tmp_h += w;
// 						 tmp_h -= shifter;
// 						 int val = (tmp_h << f) + tmp_g;
// 						 mapTable.getID(val, IDs);
// 				}
// 			}
// 		}
// 	}
// 		map<unsigned int, int>::iterator it = index.find(k1);
// 		if (it != index.end())
// 		{
// 			int tag = it->second;
// 			link_node* node = buffer[tag];
// 			if(node->weight!=0)
// 				mapTable.getID(k1, IDs); // self-loop
// 			node = node->next;
// 			while (node != NULL)
// 			{
// 				mapTable.getID(node->key, IDs);
// 				node=node->next;
// 			}
// 		}
// 		delete []tmp1;
// 		return;
// }
// void HGSS::nodePrecursorQuery(string s1, vector<string>&IDs) // same as successor query
// {
// 	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
// 	int tmp=pow(2,f)-1;
// 	unsigned short g1=hash1 & tmp;
// 	unsigned int h1 = (hash1>>f)%w;
// 	if(g1==0) g1+=1;
// 	int* tmp1 = new int[r];
// 	tmp1[0] = g1;
// 	unsigned int k1 = (h1 << f) + g1;
// 	for (int i = 1; i < r; i++)
// 	{
// 		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
// 	}
// 	for (int i = 0; i < r; i++)
// 	{
// 		int p1 = (h1 + tmp1[i]) % w;
// 		for (int k = 0; k < w; k++)
// 		{
// 			int pos = p1 + k*w;
// 			for (int j = 0; j < s; ++j)
// 			{
// 				if ((((value[pos].idx >> ((j << 3)+4))&((1 << 4) - 1)) == i) && (value[pos].dst[j] == g1))
// 				{
// 					     int tmp_g = value[pos].src[j];
// 						 int tmp_s = ((value[pos].idx >> (j << 3))&((1 << 4) - 1));

// 						 int shifter = tmp_g;
// 						 for (int v = 0; v < tmp_s; v++)
// 							 shifter = (shifter*timer + prime) % bigger_p;
// 						 int tmp_h = k;
// 						 while (tmp_h < shifter)
// 							 tmp_h += w;
// 						 tmp_h -= shifter;

// 						 int val = (tmp_h << f) + tmp_g;
// 						 mapTable.getID(val, IDs);
// 				}
// 			}
// 		}
// 	}

// 		map<unsigned int, int>::iterator it = index.find(k1);
// 		if(it!=index.end())
// 		{
// 			if(buffer[it->second]->weight!=0)
// 				mapTable.getID(k1, IDs);
// 		}

// 			for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
// 		{
// 			int tag = it->second;
// 			int src = it->first;
// 			link_node* node = buffer[tag];
// 			node = node->next;
// 			while (node != NULL)
// 			{
// 				if(node->key == k1)
// 				{
// 					mapTable.getID(src, IDs);
// 					break;
// 				}
// 				node = node->next;
// 			}
// 		}
// 		delete []tmp1;
// 		return;
// }
int HGSS::edgeQuery(string s1, string s2, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return the weight of the edge
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int *tmp1 = new int[r];
    int *tmp2 = new int[r];
    tmp1[0] = g1;
    tmp2[0] = g2;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
        tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
    }
    long key = g1 + g2;

    for (int i = 0; i < p; i++)
    {
        key = (key * timer + prime) % bigger_p;
        int index = key % (r * r);
        int index1 = index / r;
        int index2 = index % r;
        int p1 = (h1 + tmp1[index1]) % w;
        int p2 = (h2 + tmp2[index2]) % w;
        int pos = p1 * w + p2;
        for (int j = 0; j < s; j++)
        {

            if ((((value[pos].idx >> (j << 3)) & ((1 << 8) - 1)) == (index1 | (index2 << 4))) && (value[pos].src[j] == g1) && (value[pos].dst[j] == g2))
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                for (int t = 0; t < win_k; t++)
                    weight += value[pos].C[j][t];
                return weight;
            }
        }
    }
    unsigned int k1 = (h1 << f) + g1;
    unsigned int k2 = (h2 << f) + g2;
    map<unsigned int, int>::iterator it = index.find(k1);
    if (it != index.end())
    {
        int tag = it->second;
        link_node *node = buffer[tag];
        while (node != NULL)
        {
            if (node->key == k2)
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
                for (int t = b_i1; t < win_k; t++)
                    weight += node->C[t];
                return weight;
            }
            node = node->next;
        }
    }
    delete[] tmp1;
    delete[] tmp2;
    return 0;
}

int HGSS::edgeQueryWithLabel(string s1, string s2, string l_e, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return the weight of the edge
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int *tmp1 = new int[r];
    int *tmp2 = new int[r];
    tmp1[0] = g1;
    tmp2[0] = g2;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
        tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
    }
    long key = g1 + g2;

    // 引入质数分组
    int prime_group_num = 5; //一个组包含多少个质数
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    // 第一次出现，随机
    // 第二次及以后，使用第一次的随机数（否则同一个边label会对应不同的质数）
    map<string, int>::iterator it0 = random_index.find(l_e);
    // 存在label的随机数
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    // 最后一组可能会越界
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    for (int i = 0; i < p; i++)
    {
        key = (key * timer + prime) % bigger_p;
        int index = key % (r * r);
        int index1 = index / r;
        int index2 = index % r;
        int p1 = (h1 + tmp1[index1]) % w;
        int p2 = (h2 + tmp2[index2]) % w;
        int pos = p1 * w + p2;
        for (int j = 0; j < s; j++)
        {

            if ((((value[pos].idx >> (j << 3)) & ((1 << 8) - 1)) == (index1 | (index2 << 4))) && (value[pos].src[j] == g1) && (value[pos].dst[j] == g2))
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                for (int t = 0; t < win_k; t++)
                {
                    int temp_P = value[pos].P[j][t];
                    while (temp_P != 0 && temp_P % edge_label == 0)
                    {
                        weight++;
                        temp_P = temp_P / edge_label;
                    }
                }
                return weight;
            }
        }
    }
    unsigned int k1 = (h1 << f) + g1;
    unsigned int k2 = (h2 << f) + g2;
    map<unsigned int, int>::iterator it = index.find(k1);
    if (it != index.end())
    {
        int tag = it->second;
        link_node *node = buffer[tag];
        while (node != NULL)
        {
            if (node->key == k2)
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
                for (int t = b_i1; t < win_k; t++)
                {
                    int temp_P = node->P[t];
                    while (temp_P != 0 && temp_P % edge_label == 0)
                    {
                        weight++;
                        temp_P = temp_P / edge_label;
                    }
                }
                return weight;
            }
            node = node->next;
        }
    }
    delete[] tmp1;
    delete[] tmp2;
    return 0;
}

int HGSS::edgeQuery(string s1, string l_A, string s2, string l_B, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return the weight of the edge
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int *tmp1 = new int[r];
    int *tmp2 = new int[r];
    tmp1[0] = g1;
    tmp2[0] = g2;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
        tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
    }
    long key = g1 + g2;

    // 矩阵块
    unsigned int label1 = (*hfunc[0])((unsigned char *)(l_A.c_str()), l_A.length());
    unsigned int label2 = (*hfunc[0])((unsigned char *)(l_B.c_str()), l_B.length());
    unsigned int m = (label1 % (w / b));
    unsigned int n = (label2 % (w / b));

    for (int i = 0; i < p; i++)
    {
        key = (key * timer + prime) % bigger_p;
        int index = key % (r * r);
        int index1 = index / r;
        int index2 = index % r;
        int p1 = (h1 + tmp1[index1]) % b; // 在block中的位置
        int p2 = (h2 + tmp2[index2]) % b;
        int pos = (p1 + m * b) * w + (p2 + n * b);
        for (int j = 0; j < s; j++)
        {

            if ((((value[pos].idx >> (j << 3)) & ((1 << 8) - 1)) == (index1 | (index2 << 4))) && (value[pos].src[j] == g1) && (value[pos].dst[j] == g2))
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                for (int t = 0; t < win_k; t++)
                    weight += value[pos].C[j][t];
                return weight;
            }
        }
    }
    unsigned int k1 = (h1 << f) + g1;
    unsigned int k2 = (h2 << f) + g2;
    map<unsigned int, int>::iterator it = index.find(k1);
    if (it != index.end())
    {
        int tag = it->second;
        link_node *node = buffer[tag];
        while (node != NULL)
        {
            if (node->key == k2)
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
                for (int t = b_i1; t < win_k; t++)
                    weight += node->C[t];
                return weight;
            }
            node = node->next;
        }
    }
    delete[] tmp1;
    delete[] tmp2;
    return 0;
}

int HGSS::edgeQueryWithLabel(string s1, string l_A, string s2, string l_B, string l_e, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return the weight of the edge
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int *tmp1 = new int[r];
    int *tmp2 = new int[r];
    tmp1[0] = g1;
    tmp2[0] = g2;
    for (int i = 1; i < r; i++)
    {
        tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
        tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
    }
    long key = g1 + g2;

    // 矩阵块
    unsigned int label1 = (*hfunc[0])((unsigned char *)(l_A.c_str()), l_A.length());
    unsigned int label2 = (*hfunc[0])((unsigned char *)(l_B.c_str()), l_B.length());
    unsigned int m = (label1 % (w / b));
    unsigned int n = (label2 % (w / b));

    // 引入质数分组
    int prime_group_num = 5; //一个组包含多少个质数
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    // 第一次出现，随机
    // 第二次及以后，使用第一次的随机数（否则同一个边label会对应不同的质数）
    map<string, int>::iterator it0 = random_index.find(l_e);
    // 存在label的随机数
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    // 最后一组可能会越界
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    for (int i = 0; i < p; i++)
    {
        key = (key * timer + prime) % bigger_p;
        int index = key % (r * r);
        int index1 = index / r;
        int index2 = index % r;
        int p1 = (h1 + tmp1[index1]) % b; // 在block中的位置
        int p2 = (h2 + tmp2[index2]) % b;
        int pos = (p1 + m * b) * w + (p2 + n * b);
        for (int j = 0; j < s; j++)
        {
            if ((((value[pos].idx >> (j << 3)) & ((1 << 8) - 1)) == (index1 | (index2 << 4))) && (value[pos].src[j] == g1) && (value[pos].dst[j] == g2))
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                for (int t = 0; t < win_k; t++)
                {
                    int temp_P = value[pos].P[j][t];
                    while (temp_P != 0 && temp_P % edge_label == 0)
                    {
                        weight++;
                        temp_P = temp_P / edge_label;
                    }
                }
                return weight;
            }
        }
    }
    unsigned int k1 = (h1 << f) + g1;
    unsigned int k2 = (h2 << f) + g2;
    map<unsigned int, int>::iterator it = index.find(k1);
    if (it != index.end())
    {
        int tag = it->second;
        link_node *node = buffer[tag];
        while (node != NULL)
        {
            if (node->key == k2)
            {
                delete[] tmp1;
                delete[] tmp2;
                // my
                int weight = 0;
                int b_i1 = floor(difftime(end_time, node->last_T) / sub_size); // 满足时间的窗口的起始index
                for (int t = b_i1; t < win_k; t++)
                {
                    int temp_P = node->P[t];
                    while (temp_P != 0 && temp_P % edge_label == 0)
                    {
                        weight++;
                        temp_P = temp_P / edge_label;
                    }
                }
                return weight;
            }
            node = node->next;
        }
    }
    delete[] tmp1;
    delete[] tmp2;
    return 0;
}
bool HGSS::query(string s1, string s2, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return whether reachable.
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int pos;
    map<unsigned int, bool> checked;
    queue<mnode> q;
    mnode e;
    e.h = h1;
    e.g = g1;
    q.push(e);
    checked[((h1 << f) + g1)] = true;
    map<unsigned int, bool>::iterator IT;

    while (!q.empty())
    {
        e = q.front();
        h1 = e.h;
        g1 = e.g;
        int *tmp1 = new int[r];
        int *tmp2 = new int[r];
        tmp2[0] = g2;
        tmp1[0] = g1;
        for (int i = 1; i < r; i++)
        {
            tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
            tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
        }
        // 直接查找
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % w;
            for (int i2 = 0; i2 < r; i2++)
            {
                int p2 = (h2 + tmp2[i2]) % w;
                int pos = p1 * w + p2;
                for (int i3 = 0; i3 < s; i3++)
                {

                    if ((((value[pos].idx >> (i3 << 3)) & ((1 << 8) - 1)) == (i1 | (i2 << 4))) && (value[pos].src[i3] == g1) && (value[pos].dst[i3] == g2))
                    {
                        delete[] tmp1;
                        delete[] tmp2;
                        return true;
                    }
                }
            }
        }
        // 将所有出度的顶点加入queue  先buffer再matrix
        unsigned int k1 = (h1 << f) + g1;

        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            while (node != NULL)
            {
                if (node->key != k1)
                {
                    // 有效窗口内有边权重才可以当做存在，加进queue
                    bool ifvalid = false;
                    int b_i = floor(difftime(end_time, node->last_T) / sub_size);
                    for (int t = b_i; t < win_k; t++)
                    {
                        if (node->C[t] > 0)
                        {
                            ifvalid = true;
                            break;
                        }
                    }

                    if (ifvalid)
                    {
                        unsigned int val = node->key;
                        unsigned int temp_h = (val) >> f;
                        unsigned int tmp = pow(2, f);
                        unsigned short temp_g = (val % tmp);
                        if ((temp_h == h2) && (temp_g == g2))
                        {
                            delete[] tmp1;
                            delete[] tmp2;
                            return true;
                        }

                        IT = checked.find(val);
                        if (IT == checked.end())
                        {
                            mnode temp_e;
                            temp_e.h = temp_h;
                            temp_e.g = temp_g;
                            q.push(temp_e);
                            checked[val] = true;
                        }
                    }
                }
                node = node->next;
            }
        }
        //	 }
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % w;
            for (int i2 = 0; i2 < w; i2++)
            {
                int pos = p1 * w + i2;
                for (int i3 = 0; i3 < s; i3++)
                {

                    if (value[pos].src[i3] == g1 && (((value[pos].idx >> (i3 << 3)) & ((1 << 4) - 1)) == i1))
                    {
                        int tmp_g = value[pos].dst[i3];
                        int tmp_s = ((value[pos].idx >> ((i3 << 3) + 4)) & ((1 << 4) - 1));

                        int shifter = tmp_g;
                        for (int v = 0; v < tmp_s; v++)
                            shifter = (shifter * timer + prime) % bigger_p;
                        int tmp_h = i2;
                        while (tmp_h < shifter)
                            tmp_h += w;
                        tmp_h -= shifter;

                        unsigned int val = (tmp_h << f) + tmp_g;

                        IT = checked.find(val);
                        if (IT == checked.end())
                        {

                            mnode tmp_e;
                            tmp_e.h = tmp_h;
                            tmp_e.g = tmp_g;
                            q.push(tmp_e);
                            checked[val] = true;
                        }
                    }
                }
            }
        }
        delete[] tmp1;
        delete[] tmp2;
        q.pop();
    }
    return false;
}

bool HGSS::queryWithLabel(string s1, string s2, string l_e, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return whether reachable.
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int pos;
    map<unsigned int, bool> checked;
    queue<mnode> q;
    mnode e;
    e.h = h1;
    e.g = g1;
    q.push(e);
    checked[((h1 << f) + g1)] = true;
    map<unsigned int, bool>::iterator IT;

    // 引入质数分组
    int prime_group_num = 5; //一个组包含多少个质数
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    // 第一次出现，随机
    // 第二次及以后，使用第一次的随机数（否则同一个边label会对应不同的质数）
    map<string, int>::iterator it0 = random_index.find(l_e);
    // 存在label的随机数
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    // 最后一组可能会越界
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    while (!q.empty())
    {
        e = q.front();
        h1 = e.h;
        g1 = e.g;
        int *tmp1 = new int[r];
        int *tmp2 = new int[r];
        tmp2[0] = g2;
        tmp1[0] = g1;
        for (int i = 1; i < r; i++)
        {
            tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
            tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
        }
        // 直接查找
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % w;
            for (int i2 = 0; i2 < r; i2++)
            {
                int p2 = (h2 + tmp2[i2]) % w;
                int pos = p1 * w + p2;
                for (int i3 = 0; i3 < s; i3++)
                {

                    if ((((value[pos].idx >> (i3 << 3)) & ((1 << 8) - 1)) == (i1 | (i2 << 4))) && (value[pos].src[i3] == g1) && (value[pos].dst[i3] == g2))
                    {
                        delete[] tmp1;
                        delete[] tmp2;
                        // 找是否有对应的边label
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[i3][t];
                            if (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                return true;
                            }
                        }
                        // return true;
                    }
                }
            }
        }
        // 将所有出度的顶点加入queue  先buffer再matrix
        unsigned int k1 = (h1 << f) + g1;

        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            while (node != NULL)
            {
                if (node->key != k1)
                {
                    // 有效窗口内有边权重才可以当做存在，加进queue
                    bool ifvalid = false;
                    int b_i = floor(difftime(end_time, node->last_T) / sub_size);
                    for (int t = b_i; t < win_k; t++)
                    {
                        if (node->P[t] % edge_label == 0) //label限制
                        {
                            ifvalid = true;
                            break;
                        }
                    }

                    if (ifvalid)
                    {
                        unsigned int val = node->key;
                        unsigned int temp_h = (val) >> f;
                        unsigned int tmp = pow(2, f);
                        unsigned short temp_g = (val % tmp);
                        if ((temp_h == h2) && (temp_g == g2))
                        {
                            delete[] tmp1;
                            delete[] tmp2;
                            return true;
                        }

                        IT = checked.find(val);
                        if (IT == checked.end())
                        {
                            mnode temp_e;
                            temp_e.h = temp_h;
                            temp_e.g = temp_g;
                            q.push(temp_e);
                            checked[val] = true;
                        }
                    }
                }
                node = node->next;
            }
        }
        //	 }
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % w;
            for (int i2 = 0; i2 < w; i2++)
            {
                int pos = p1 * w + i2;
                for (int i3 = 0; i3 < s; i3++)
                {

                    if (value[pos].src[i3] == g1 && (((value[pos].idx >> (i3 << 3)) & ((1 << 4) - 1)) == i1))
                    {
                        int tmp_g = value[pos].dst[i3];
                        int tmp_s = ((value[pos].idx >> ((i3 << 3) + 4)) & ((1 << 4) - 1));

                        int shifter = tmp_g;
                        for (int v = 0; v < tmp_s; v++)
                            shifter = (shifter * timer + prime) % bigger_p;
                        int tmp_h = i2;
                        while (tmp_h < shifter)
                            tmp_h += w;
                        tmp_h -= shifter;

                        // label限制
                        bool f = false;
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[i3][t];
                            if (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                f = true;
                                break;
                            }
                        }

                        if (f)
                        {
                            unsigned int val = (tmp_h << f) + tmp_g;

                            IT = checked.find(val);
                            if (IT == checked.end())
                            {

                                mnode tmp_e;
                                tmp_e.h = tmp_h;
                                tmp_e.g = tmp_g;
                                q.push(tmp_e);
                                checked[val] = true;
                            }
                        }
                    }
                }
            }
        }
        delete[] tmp1;
        delete[] tmp2;
        q.pop();
    }
    return false;
}

bool HGSS::query(string s1, string l_A, string s2, string l_B, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return whether reachable.
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int pos;
    map<unsigned int, bool> checked;
    queue<mnode> q;
    mnode e;
    e.h = h1;
    e.g = g1;
    q.push(e);
    checked[((h1 << f) + g1)] = true;
    map<unsigned int, bool>::iterator IT;

    // 矩阵块
    unsigned int label1 = (*hfunc[0])((unsigned char *)(l_A.c_str()), l_A.length());
    unsigned int label2 = (*hfunc[0])((unsigned char *)(l_B.c_str()), l_B.length());
    unsigned int m = (label1 % (w / b));
    unsigned int n = (label2 % (w / b));

    while (!q.empty())
    {
        e = q.front();
        h1 = e.h;
        g1 = e.g;
        unsigned int k1 = (h1 << f) + g1;
        // 更新m
        string label_tmp = k2label[k1];
        label1 = (*hfunc[0])((unsigned char *)(label_tmp.c_str()), label_tmp.length());
        m = (label1 % (w / b));

        int *tmp1 = new int[r];
        int *tmp2 = new int[r];
        tmp2[0] = g2;
        tmp1[0] = g1;
        for (int i = 1; i < r; i++)
        {
            tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
            tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
        }
        // 直接查找
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % b; // block
            for (int i2 = 0; i2 < r; i2++)
            {
                int p2 = (h2 + tmp2[i2]) % b; // block
                int pos = (p1 + m * b) * w + (p2 + n * b);
                for (int i3 = 0; i3 < s; i3++)
                {

                    if ((((value[pos].idx >> (i3 << 3)) & ((1 << 8) - 1)) == (i1 | (i2 << 4))) && (value[pos].src[i3] == g1) && (value[pos].dst[i3] == g2))
                    {
                        delete[] tmp1;
                        delete[] tmp2;
                        return true;
                    }
                }
            }
        }

        // 将所有出度的顶点加入queue  先buffer再matrix
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            while (node != NULL)
            {
                if (node->key != k1)
                {
                    // 有效窗口内有边权重才可以当做存在，加进queue
                    bool ifvalid = false;
                    int b_i = floor(difftime(end_time, node->last_T) / sub_size);
                    for (int t = b_i; t < win_k; t++)
                    {
                        if (node->C[t] > 0)
                        {
                            ifvalid = true;
                            break;
                        }
                    }

                    if (ifvalid)
                    {
                        unsigned int val = node->key;
                        unsigned int temp_h = (val) >> f;
                        unsigned int tmp = pow(2, f);
                        unsigned short temp_g = (val % tmp);
                        if ((temp_h == h2) && (temp_g == g2))
                        {
                            delete[] tmp1;
                            delete[] tmp2;
                            return true;
                        }

                        IT = checked.find(val);
                        if (IT == checked.end())
                        {
                            mnode temp_e;
                            temp_e.h = temp_h;
                            temp_e.g = temp_g;
                            q.push(temp_e);
                            checked[val] = true;
                        }
                    }
                }
                node = node->next;
            }
        }
        //	 }
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % b; // block
            for (int i2 = 0; i2 < w; i2++)
            {
                int pos = (m * b + p1) * w + i2;
                for (int i3 = 0; i3 < s; i3++)
                {

                    if (value[pos].src[i3] == g1 && (((value[pos].idx >> (i3 << 3)) & ((1 << 4) - 1)) == i1))
                    {
                        int tmp_g = value[pos].dst[i3];
                        // 取index2
                        int tmp_s = ((value[pos].idx >> ((i3 << 3) + 4)) & ((1 << 4) - 1));
                        int tmp_h = value[pos].h2[i3];

                        unsigned int val = (tmp_h << f) + tmp_g;
                        IT = checked.find(val);
                        if (IT == checked.end())
                        {

                            mnode tmp_e;
                            tmp_e.h = tmp_h;
                            tmp_e.g = tmp_g;
                            q.push(tmp_e);
                            checked[val] = true;
                        }
                    }
                }
            }
        }
        delete[] tmp1;
        delete[] tmp2;
        q.pop();
    }
    return false;
}

bool HGSS::queryWithLabel(string s1, string l_A, string s2, string l_B, string l_e, time_t end_time) // s1 is the ID of the source node, s2 is the ID of the destination node, return whether reachable.
{
    unsigned int hash1 = (*hfunc[0])((unsigned char *)(s1.c_str()), s1.length());
    unsigned int hash2 = (*hfunc[0])((unsigned char *)(s2.c_str()), s2.length());
    int tmp = pow(2, f) - 1;
    unsigned short g1 = hash1 & tmp;
    if (g1 == 0)
        g1 += 1;
    unsigned int h1 = (hash1 >> f) % w;
    unsigned short g2 = hash2 & tmp;
    if (g2 == 0)
        g2 += 1;
    unsigned int h2 = (hash2 >> f) % w;
    int pos;
    map<unsigned int, bool> checked;
    queue<mnode> q;
    mnode e;
    e.h = h1;
    e.g = g1;
    q.push(e);
    checked[((h1 << f) + g1)] = true;
    map<unsigned int, bool>::iterator IT;

    // 矩阵块
    unsigned int label1 = (*hfunc[0])((unsigned char *)(l_A.c_str()), l_A.length());
    unsigned int label2 = (*hfunc[0])((unsigned char *)(l_B.c_str()), l_B.length());
    unsigned int m = (label1 % (w / b));
    unsigned int n = (label2 % (w / b));

    // 引入质数分组
    int prime_group_num = 5; //一个组包含多少个质数
    int group_s = ceil((float)prime_num / (float)prime_group_num);
    unsigned int labelhash = (*hfunc[0])((unsigned char *)(l_e.c_str()), l_e.length());
    int prime_group = labelhash % group_s;
    int r_number;
    // 第一次出现，随机
    // 第二次及以后，使用第一次的随机数（否则同一个边label会对应不同的质数）
    map<string, int>::iterator it0 = random_index.find(l_e);
    // 存在label的随机数
    if (it0 != random_index.end())
        r_number = it0->second;
    else
    {
        srand(labelhash);
        r_number = (rand() % prime_group_num);
        random_index[l_e] = r_number;
    }
    unsigned int edge_label_index = prime_group * prime_group_num + r_number;
    // 最后一组可能会越界
    if (edge_label_index > prime_num)
        edge_label_index = prime_num - 1;
    int edge_label = Prime[edge_label_index];

    while (!q.empty())
    {
        e = q.front();
        h1 = e.h;
        g1 = e.g;
        unsigned int k1 = (h1 << f) + g1;
        // 更新m
        string label_tmp = k2label[k1];
        label1 = (*hfunc[0])((unsigned char *)(label_tmp.c_str()), label_tmp.length());
        m = (label1 % (w / b));

        int *tmp1 = new int[r];
        int *tmp2 = new int[r];
        tmp2[0] = g2;
        tmp1[0] = g1;
        for (int i = 1; i < r; i++)
        {
            tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
            tmp2[i] = (tmp2[i - 1] * timer + prime) % bigger_p;
        }
        // 直接查找
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % b; // block
            for (int i2 = 0; i2 < r; i2++)
            {
                int p2 = (h2 + tmp2[i2]) % b; // block
                int pos = (p1 + m * b) * w + (p2 + n * b);
                for (int i3 = 0; i3 < s; i3++)
                {

                    if ((((value[pos].idx >> (i3 << 3)) & ((1 << 8) - 1)) == (i1 | (i2 << 4))) && (value[pos].src[i3] == g1) && (value[pos].dst[i3] == g2))
                    {

                        // 找是否有对应的边label
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[i3][t];
                            if (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                delete[] tmp1;
                                delete[] tmp2;
                                return true;
                            }
                        }
                        // return true;
                    }
                }
            }
        }

        // 将所有出度的顶点加入queue  先buffer再matrix
        map<unsigned int, int>::iterator it = index.find(k1);
        if (it != index.end())
        {
            int tag = it->second;
            link_node *node = buffer[tag];
            while (node != NULL)
            {
                if (node->key != k1)
                {
                    // 有效窗口内有边权重才可以当做存在，加进queue
                    bool ifvalid = false;
                    int b_i = floor(difftime(end_time, node->last_T) / sub_size);
                    for (int t = b_i; t < win_k; t++)
                    {
                        if (node->P[t] % edge_label == 0) //label限制
                        {
                            ifvalid = true;
                            break;
                        }
                    }

                    if (ifvalid)
                    {
                        unsigned int val = node->key;
                        unsigned int temp_h = (val) >> f;
                        unsigned int tmp = pow(2, f);
                        unsigned short temp_g = (val % tmp);
                        if ((temp_h == h2) && (temp_g == g2))
                        {
                            delete[] tmp1;
                            delete[] tmp2;
                            return true;
                        }

                        IT = checked.find(val);
                        if (IT == checked.end())
                        {
                            mnode temp_e;
                            temp_e.h = temp_h;
                            temp_e.g = temp_g;
                            q.push(temp_e);
                            checked[val] = true;
                        }
                    }
                }
                node = node->next;
            }
        }
        //	 }
        for (int i1 = 0; i1 < r; i1++)
        {
            int p1 = (h1 + tmp1[i1]) % b; // block
            for (int i2 = 0; i2 < w; i2++)
            {
                int pos = (m * b + p1) * w + i2;
                for (int i3 = 0; i3 < s; i3++)
                {

                    if (value[pos].src[i3] == g1 && (((value[pos].idx >> (i3 << 3)) & ((1 << 4) - 1)) == i1))
                    {
                        int tmp_g = value[pos].dst[i3];
                        // 取index2
                        int tmp_s = ((value[pos].idx >> ((i3 << 3) + 4)) & ((1 << 4) - 1));
                        int tmp_h = value[pos].h2[i3];

                        // label限制
                        bool f = false;
                        for (int t = 0; t < win_k; t++)
                        {
                            int temp_P = value[pos].P[i3][t];
                            if (temp_P != 0 && temp_P % edge_label == 0)
                            {
                                f = true;
                                break;
                            }
                        }
                        if (f)
                        {
                            unsigned int val = (tmp_h << f) + tmp_g;
                            IT = checked.find(val);
                            if (IT == checked.end())
                            {

                                mnode tmp_e;
                                tmp_e.h = tmp_h;
                                tmp_e.g = tmp_g;
                                q.push(tmp_e);
                                checked[val] = true;
                            }
                        }
                    }
                }
            }
        }
        delete[] tmp1;
        delete[] tmp2;
        q.pop();
    }
    return false;
}

// /*type 0 is for successor query, type 1 is for precusor query*/
// int HGSS::nodeDegreeQuery(string s1, int type) // s1 is the ID of the queried node, return the in/out degree
// {
// 	int degree = 0;
// 	unsigned int hash1 = (*hfunc[0])((unsigned char*)(s1.c_str()), s1.length());
// 	int tmp = pow(2, f) - 1;
// 	unsigned short g1 = hash1 & tmp;
// 	if (g1 == 0) g1 += 1;
// 	unsigned int h1 = (hash1 >> f) % w;
// 	int* tmp1 = new int[r];
// 	tmp1[0] = g1;
// 	for (int i = 1; i < r; i++)
// 	{
// 		tmp1[i] = (tmp1[i - 1] * timer + prime) % bigger_p;
// 	}
// 	for (int i = 0; i < r; i++)
// 	{
// 		int p1 = (h1 + tmp1[i]) % w;
// 		for (int k = 0; k < w; k++)
// 		{
// 			if (type == 0)/*successor query*/
// 			{
// 				int pos = p1*w + k;
// 				for (int j = 0; j < s; ++j)
// 				{
// 					if ((((value[pos].idx >> ((j << 3)))&((1 << 4) - 1)) == i) && (value[pos].src[j] == g1))
// 					{
// 						 int tmp_g = value[pos].dst[j];
// 						 int tmp_s = ((value[pos].idx >> ((j << 3) + 4))&((1 << 4) - 1));

// 						 int shifter = tmp_g;
// 						 for (int v = 0; v < tmp_s; v++)
// 							 shifter = (shifter*timer + prime) % bigger_p;
// 						 int tmp_h = k;
// 						 while (tmp_h < shifter)
// 							 tmp_h += w;
// 						 tmp_h -= shifter;
// 						 unsigned int val = (tmp_h << f) + tmp_g;

// 						 degree+=mapTable.countIDnums(val);
// 					}
// 				}
// 			}
// 			else if (type == 1)/*precursor query*/
// 			{
// 				int pos = p1 + k*w;
// 				for (int j = 0; j < s; ++j)
// 				{
// 					if ((((value[pos].idx >> ((j << 3) + 4))&((1 << 4) - 1)) == i) && (value[pos].dst[j] == g1))
// 					{
// 						 int tmp_g = value[pos].src[j];
// 						 int tmp_s = ((value[pos].idx >> (j << 3))&((1 << 4) - 1));

// 						 int shifter = tmp_g;
// 						 for (int v = 0; v < tmp_s; v++)
// 							 shifter = (shifter*timer + prime) % bigger_p;
// 						 int tmp_h = k;
// 						 while (tmp_h < shifter)
// 							 tmp_h += w;
// 						 tmp_h -= shifter;
// 						 unsigned int val = (tmp_h << f) + tmp_g;

// 						degree+=mapTable.countIDnums(val);
// 					}
// 				}
// 			}
// 		}
// 	}

// 	if (type == 0)
// 	{
// 		unsigned int k1 = (h1 << f) + g1;
// 		map<unsigned int, int>::iterator it = index.find(k1);
// 		if (it != index.end())
// 		{
// 			int tag = it->second;
// 			link_node* node = buffer[tag];
// 			if(node->weight!=0);
// 				degree += mapTable.countIDnums(k1); // address self-loops first
// 			node = node->next;
// 			while (node != NULL)
// 			{
// 				degree+=mapTable.countIDnums(node->key);
// 				node = node->next;
// 			}
// 		}
// 	}
// 	else if (type == 1)
// 	{
// 		unsigned int k1 = (h1 << f) + g1;
// 		map<unsigned int, int>::iterator it = index.find(k1); // address self-loops first
// 		if(it!=index.end())
// 		{
// 			if(buffer[it->second]->weight!=0)
// 				degree += mapTable.countIDnums(k1);
// 		}
// 		for (map<unsigned int, int>::iterator it = index.begin(); it != index.end(); ++it)
// 		{
// 			int tag = it->second;
// 			link_node* node = buffer[tag];
// 			unsigned int src=node->key;
// 			node = node->next;
// 			while (node != NULL)
// 			{
// 				if (node->key == k1)
// 				{
// 					degree+=mapTable.countIDnums(src);
// 					break;
// 				 }
// 				node = node->next;
// 			}
// 		}
// 	}
// 	delete[]tmp1;
// 	return degree;
// }

// void hvinsert(unsigned int hash, std::vector<hashvalue> &v)
// {
// 	bool find=false;
// 	for(int i=0;i<v.size();i++)
// 	{
// 		if (hash==v[i].key)
// 		{
// 			v[i].IDnum++;
// 			find=true;
// 			break;
// 		}
// 	}
// 	if(!find)
// 	{
// 		hashvalue hv;
// 		hv.key=hash;
// 		hv.IDnum=1;
// 		v.push_back(hv);
// 	}
// 	return;
// }

// int HGSS::TriangleCounting()
// {
// 	GSketch* gs = new GSketch;
// 	for(int i=0;i<tablesize;i++)
// 	{
// 		hashTableNode<string> *np=mapTable.table[i];
// 		vector<hashvalue> nodes;
// 		for(;np!=NULL;np=np->next)
// 		 	hvinsert(np->key,nodes);
// 		int nodenum=nodes.size();
// 		for(int j=0;j<nodenum;j++)
// 		{
// 			//cout<<nodes[j].key<<' '<<nodes[j].IDnum<<endl;
// 			gs->insert_node(nodes[j].key, nodes[j].IDnum);
// 		}
// 		nodes.clear();
// 	}
// 	for(int i=0;i<w;i++)
// 	{
// 		int row=i*w;
// 		for(int k=0;k<w;k++)
// 		{
// 			int pos=row+k;
// 			for(int j=0;j<s;j++)
// 			{
// 				if(value[pos].src[j]>0)
// 				{
// 					unsigned int tmp_g1 = value[pos].src[j];
// 					unsigned int tmp_s1 = ((value[pos].idx >> (j << 3))&((1 << 4) - 1));

// 					int shifter = tmp_g1;
// 					for (int v = 0; v < tmp_s1; v++)
// 						shifter = (shifter*timer + prime) % bigger_p;
// 					int tmp_h1 = i;
// 					while (tmp_h1 < shifter)
// 						tmp_h1 += w;
// 					tmp_h1 -= shifter;

// 					unsigned int val1 = (tmp_h1 << f) + tmp_g1;

// 					 int tmp_g2 = value[pos].dst[j];
// 					 int tmp_s2 = ((value[pos].idx >> ((j << 3) + 4))&((1 << 4) - 1));

// 					 shifter = tmp_g2;
// 					 for (int v = 0; v < tmp_s2; v++)
// 						shifter = (shifter*timer + prime) % bigger_p;
// 					  int tmp_h2 = k;
// 					  while (tmp_h2 < shifter)
// 							 tmp_h2 += w;
// 					  tmp_h2 -= shifter;
// 					  unsigned int val2 = (tmp_h2 << f) + tmp_g2;
// 					  gs->insert_edge(val1,val2);
// 				}
// 			}
// 		}
// 	}
// 	for(int i=0;i<buffer.size();i++)
// 	{
// 		link_node* np=buffer[i];
// 		unsigned int src=np->key;
// 		np=np->next;
// 		for(;np!=NULL;np=np->next)
// 		{
// 		  unsigned int dst=np->key;
// 		  gs->insert_edge(src, dst);
// 		}
// 	}
// 	gs->clean();
// 	return gs->countTriangle();
// 	delete gs;
// }
