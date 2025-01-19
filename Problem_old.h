#ifndef PROBLEM_H
#define PROBLEM_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <unordered_map>
#include <chrono>
#include <set>
#include <algorithm>
#include <unordered_set>
#include "Graph.h"

using namespace std;
using namespace chrono;

struct HeuristicaConstrutivaResponse
{
    int num_cores_encontradas;
    int num_vertices;
    int num_grau_maximo;
    int delta;
    int num_arestas;
    int densidade_arestas;
    int size_clique;
    float delta_clique;
};

class Problem
{
public:
    int calcularConflitos(Graph &graph, vector<int> &cores);
    int calcularConflitos(Graph &graph, vector<int> &cores, int vertice, int cor);
    bool tabucol(Graph &graph, int max_iter, int n_cores);
    bool tabucol2(Graph &graph, int max_iter, int n_cores);
    bool tabucol3(Graph &graph, int max_iter, int n_cores);
    bool buscaClique(Graph &graph, int k, vector<int> &clique, int comeco);
    bool temKClique(Graph &graph, int k);
    bool temKCliqueBK(Graph &graph, int k);
    void bronKerbosch(Graph &graph, vector<int> &clique,
                      unordered_set<int> &candidatos, unordered_set<int> &excluidos, int k, bool &achou);
    vector<int> greedyClique(Graph &graph);
    HeuristicaConstrutivaResponse heuristicaConstrutiva(Graph &graph, int k);
};

vector<int> Problem::greedyClique(Graph &graph)
{
    int n = graph.get_num_vertices();
    vector<int> vertices(n);
    for (int i = 0; i < n; ++i)
        vertices[i] = i;

    // Ordena os vértices pelo grau decrescente
    sort(vertices.begin(), vertices.end(), [&](int a, int b)
         { return graph.get_neighbors(a).size() > graph.get_neighbors(b).size(); });

    vector<int> clique;
    for (int v : vertices)
    {
        bool canAdd = true;
        for (int u : clique)
        {
            if (!graph.ehVizinho(u, v))
            {
                canAdd = false;
                break;
            }
        }
        if (canAdd)
            clique.push_back(v);
    }
    return clique;
}

void Problem::bronKerbosch(Graph &graph, vector<int> &clique,
                           unordered_set<int> &candidatos, unordered_set<int> &excluidos, int k, bool &achou)
{
    if (achou)
        return;

    if (candidatos.empty() && excluidos.empty())
    {
        if (clique.size() == k)
        {
            achou = true; // Encontrou um clique de tamanho k
        }
        return;
    }

    int pivo = -1, maxIntersecao = -1;
    for (int u : candidatos)
    {
        int intersecao = 0;
        for (int v : candidatos)
        {
            if (graph.ehVizinho(u, v))
            {
                ++intersecao;
            }
        }
        if (intersecao > maxIntersecao)
        {
            pivo = u;
            maxIntersecao = intersecao;
        }
    }

    unordered_set<int> candidatosCopia = candidatos;
    for (int v : candidatosCopia)
    {
        if (pivo != -1 && graph.ehVizinho(pivo, v))
        {
            continue;
        }
        // Adiciona v ao clique atual
        clique.push_back(v);

        // Calcula novos candidatos e excluídos
        unordered_set<int> newcandidatos, newexcluidos;
        for (int u : candidatos)
        {
            if (graph.ehVizinho(v, u))
            {
                newcandidatos.insert(u);
            }
        }
        for (int u : excluidos)
        {
            if (graph.ehVizinho(v, u))
            {
                newexcluidos.insert(u);
            }
        }

        // Chamada recursiva
        bronKerbosch(graph, clique, newcandidatos, newexcluidos, k, achou);

        // Remove v do clique (backtracking)
        clique.pop_back();

        // Move v do conjunto de candidatos para o conjunto excluído
        candidatos.erase(v);
        excluidos.insert(v);

        if (achou)
            return; // Interrompe se o clique foi encontrado
    }
};

bool Problem::temKCliqueBK(Graph &graph, int k)
{
    unordered_set<int> candidatos, excluidos;
    vector<int> clique;
    bool achou = false;

    for (int i = 0; i < graph.get_num_vertices(); ++i)
    {
        candidatos.insert(i);
    }

    bronKerbosch(graph, clique, candidatos, excluidos, k, achou);

    return achou;
}

bool Problem::buscaClique(Graph &graph, int k, vector<int> &clique, int comeco)
{
    cout << "clique: ";
    for (int i = 0; i < clique.size(); i++)
    {
        cout << clique[i] << " ";
    }
    cout << endl;

    if (clique.size() == k)
    {
        int n = clique.size();
        for (int i = 0; i < n; ++i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                bool achou1 = false;
                bool achou2 = false;

                for (CoupleVertice vizinho : graph.get_neighbors(clique[i]))
                {
                    if (vizinho.vertice2 == clique[j])
                    {
                        achou1 = true;
                        break;
                    }
                }

                if (!achou1)
                {
                    return false;
                }

                for (CoupleVertice vizinho : graph.get_neighbors(clique[j]))
                {
                    if (vizinho.vertice2 == clique[i])
                    {
                        achou2 = true;
                        break;
                    }
                }

                if (!achou2)
                {
                    return false;
                }
            }
        }
        return true;
    }

    for (int i = comeco; i < graph.get_num_vertices(); ++i)
    {
        clique.push_back(i);
        if (buscaClique(graph, k, clique, i + 1))
        {
            return true;
        }
        clique.pop_back();
    }
    return false;
}

bool Problem::temKClique(Graph &graph, int k)
{
    vector<int> clique;
    return buscaClique(graph, k, clique, 0);
};

HeuristicaConstrutivaResponse Problem::heuristicaConstrutiva(Graph &graph, int initialK = 10)
{

    // bool temClique = false;
    // if (temKCliqueBK(graph, k))
    // {
    //     temClique = true;
    // }
    vector<int> clique = greedyClique(graph);
    int k = clique.size();

    // while (!temKCliqueBK(graph, k))
    // {
    //     k--;
    // }

    int num_arestas = 0;
    for (int i = 0; i < graph.get_num_vertices(); i++)
    {
        num_arestas += graph.get_neighbors(i).size();
    }

    int num_vertices = graph.get_num_vertices();

    vector<int> cores(num_vertices, -1);
    vector<int> saturacao(num_vertices, 0);
    vector<int> graus(num_vertices, 0);

    int grauMaximo = 0;

    for (int i = 0; i < num_vertices; i++)
    {
        graus[i] = graph.get_neighbors(i).size();
        if (graus[i] > grauMaximo)
        {
            grauMaximo = graus[i];
        }
    }

    for (int passada = 0; passada < num_vertices; passada++)
    {
        int melhorVertice = -1;
        int maiorSaturacao = -1;
        int maiorGrau = -1;

        for (int vertice = 0; vertice < num_vertices; vertice++)
        {
            if (cores[vertice] == -1)
            {
                if (saturacao[vertice] > maiorSaturacao || (saturacao[vertice] == maiorSaturacao && graus[vertice] > maiorGrau))
                {
                    melhorVertice = vertice;
                    maiorSaturacao = saturacao[vertice];
                    maiorGrau = graus[vertice];
                }
            }
        }

        int menorCor = 1;
        set<int> coresIndisponiveis;
        for (CoupleVertice vizinho : graph.get_neighbors(melhorVertice))
        {
            if (cores[vizinho.vertice2] != -1)
            {
                coresIndisponiveis.insert(cores[vizinho.vertice2]);
            }
        }

        while (coresIndisponiveis.count(menorCor) > 0)
        {
            menorCor++;
        }

        cores[melhorVertice] = menorCor;
        for (CoupleVertice vizinho : graph.get_neighbors(melhorVertice))
        {
            saturacao[vizinho.vertice2]++;
        }
    }

    int maxCor = 0;
    for (int i = 0; i < num_vertices; i++)
    {
        if (cores[i] > maxCor)
        {
            maxCor = cores[i];
        }
    }

    HeuristicaConstrutivaResponse response;

    response.size_clique = k;
    response.delta_clique = (float)(maxCor - k) / k;
    response.num_cores_encontradas = maxCor;
    response.num_grau_maximo = grauMaximo;
    response.delta = grauMaximo - maxCor;
    response.num_arestas = num_arestas;
    response.densidade_arestas = (float)(num_arestas) / num_vertices;
    response.num_vertices = num_vertices;

    return response;
}

int Problem::calcularConflitos(Graph &graph, vector<int> &cores)
{
    int conflitos = 0;
    for (int i = 0; i < graph.get_num_vertices(); i++)
    {
        for (int j = 0; j < graph.get_neighbors(i).size(); j++)
        {
            int vertice2 = graph.get_neighbors(i)[j].vertice2;
            if (cores[i] == cores[vertice2])
            {
                conflitos++;
            }
        }
    }
    return conflitos;
}

int Problem::calcularConflitos(Graph &graph, vector<int> &cores, int vertice, int cor)
{
    int conflitos = 0;
    for (int j = 0; j < graph.get_neighbors(vertice).size(); j++)
    {
        int vertice2 = graph.get_neighbors(vertice)[j].vertice2;
        if (cor == cores[vertice2])
        {
            conflitos++;
        }
    }
    return conflitos;
};

bool Problem::tabucol(Graph &graph, int max_iter, int n_cores)
{
    auto inicio = high_resolution_clock::now();
    int n_vertices = graph.get_num_vertices();
    cout << "n_vertices: " << n_vertices << endl;
    vector<int> cores(n_vertices);
    vector<vector<int>> lista_tabu(n_vertices, vector<int>(n_cores, 0));
    int iteracao = 0;
    int melhorConflito = numeric_limits<int>::max();

    srand(time(0));
    for (int i = 0; i < n_vertices; i++)
    {
        cores[i] = rand() % n_cores;
    }

    int conflitos = calcularConflitos(graph, cores);
    melhorConflito = conflitos;

    while (conflitos > 0 && iteracao < max_iter)
    {
        cout << "Iteracao: " << iteracao << " / Conflitos: " << conflitos << " / Melhor conflito: " << melhorConflito << endl;
        int melhorVizinho = -1;
        int melhorCor = -1;
        int melhorDelta = numeric_limits<int>::max();

        for (int vertice = 0; vertice < n_vertices; vertice++)
        {
            int corAtual = cores[vertice];
            for (int novaCor = 0; novaCor < n_cores; novaCor++)
            {
                if (novaCor != corAtual && lista_tabu[vertice][novaCor] <= iteracao)
                {
                    cores[vertice] = novaCor;
                    int delta = calcularConflitos(graph, cores) - conflitos;
                    cores[vertice] = corAtual;

                    if (delta < melhorDelta)
                    {
                        melhorDelta = delta;
                        melhorVizinho = vertice;
                        melhorCor = novaCor;
                    }
                }
            }
        }

        if (melhorVizinho != -1)
        {
            cores[melhorVizinho] = melhorCor;
            conflitos += melhorDelta;
            lista_tabu[melhorVizinho][melhorCor] = iteracao + n_vertices / 2;
            if (conflitos < melhorConflito)
            {
                melhorConflito = conflitos;
            }
        }
        iteracao++;
    }

    auto fim = high_resolution_clock::now();
    auto duracao = duration_cast<milliseconds>(fim - inicio).count();

    if (conflitos == 0)
    {
        cout << "Solucao otima encontrada com " << n_cores << " cores" << endl;
    }
    else
    {
        cout << "Solucao nao otima encontrada" << endl;
    }

    cout << "Tempo de execucao: " << duracao << "ms ";
    cout << "/ Iterações: " << iteracao << endl;

    return conflitos == 0;
};

bool Problem::tabucol2(Graph &graph, int max_iter, int n_cores)
{
    auto inicio = high_resolution_clock::now();
    int n_vertices = graph.get_num_vertices();
    cout << "n_vertices: " << n_vertices << endl;
    vector<int> cores(n_vertices);
    vector<vector<int>> lista_tabu(n_vertices, vector<int>(n_cores, 0));
    int iteracao = 0;
    int melhorConflito = numeric_limits<int>::max();

    srand(time(0));
    int conflitos = 0;
    for (int i = 0; i < n_vertices; i++)
    {
        cores[i] = rand() % n_cores;
        conflitos += calcularConflitos(graph, cores, i, cores[i]);
    }
    conflitos /= 2;
    melhorConflito = conflitos;

    while (conflitos > 0 && iteracao < max_iter)
    {

        int melhorVizinho = -1;
        int melhorCor = -1;
        int melhorDelta = numeric_limits<int>::max();

        for (int vertice = 0; vertice < n_vertices; vertice++)
        {
            int conflitos_vertice = calcularConflitos(graph, cores, vertice, cores[vertice]);
            if (conflitos_vertice > 0)
            {
                int corAtual = cores[vertice];
                for (int novaCor = 0; novaCor < n_cores; novaCor++)
                {
                    if (novaCor != corAtual && lista_tabu[vertice][novaCor] <= iteracao)
                    {

                        int delta = calcularConflitos(graph, cores, vertice, novaCor) - conflitos_vertice;

                        if (delta < melhorDelta)
                        {
                            melhorDelta = delta;
                            melhorVizinho = vertice;
                            melhorCor = novaCor;
                        }
                    }
                }
            }
        }

        if (melhorVizinho != -1)
        {
            int corAnterior = cores[melhorVizinho];
            cores[melhorVizinho] = melhorCor;
            conflitos += melhorDelta;
            lista_tabu[melhorVizinho][corAnterior] = iteracao + n_vertices / 2;
            if (conflitos < melhorConflito)
            {
                melhorConflito = conflitos;
            }
        }
        iteracao++;
    }

    auto fim = high_resolution_clock::now();
    auto duracao = duration_cast<milliseconds>(fim - inicio).count();

    if (conflitos <= 0)
    {
        cout << "Solucao otima encontrada com " << n_cores << " cores" << endl;
    }
    else
    {
        cout << "Solucao nao otima encontrada" << endl;
    }

    cout << "Tempo de execucao: " << duracao << "ms ";
    cout << "/ Iterações: " << iteracao << endl;

    return conflitos <= 0;
};

bool Problem::tabucol3(Graph &graph, int max_iter, int n_cores)
{
    auto inicio = high_resolution_clock::now();
    int n_vertices = graph.get_num_vertices();
    cout << "n_vertices: " << n_vertices << endl;
    vector<int> cores(n_vertices);
    vector<vector<int>> lista_tabu(n_vertices, vector<int>(n_cores, 0));
    int iteracao = 0;
    int melhorConflito = numeric_limits<int>::max();

    srand(time(0));
    int conflitos = 0;
    for (int i = 0; i < n_vertices; i++)
    {
        cores[i] = rand() % n_cores;
        conflitos += calcularConflitos(graph, cores, i, cores[i]);
    }
    conflitos /= 2;
    melhorConflito = conflitos;

    while (conflitos > 0 && iteracao < max_iter)
    {

        int melhorVizinho = -1;
        int melhorCor = -1;
        int melhorDelta = numeric_limits<int>::max();

        for (int vertice = 0; vertice < n_vertices; vertice++)
        {
            int conflitos_vertice = calcularConflitos(graph, cores, vertice, cores[vertice]);
            if (conflitos_vertice > 0)
            {
                int corAtual = cores[vertice];
                for (int novaCor = 0; novaCor < n_cores; novaCor++)
                {
                    if (novaCor != corAtual && lista_tabu[vertice][novaCor] <= iteracao)
                    {

                        int delta = calcularConflitos(graph, cores, vertice, novaCor) - conflitos_vertice;

                        if (delta < melhorDelta)
                        {
                            melhorDelta = delta;
                            melhorVizinho = vertice;
                            melhorCor = novaCor;
                        }
                    }
                }
            }
        }

        if (melhorVizinho != -1)
        {
            int corAnterior = cores[melhorVizinho];
            cores[melhorVizinho] = melhorCor;
            conflitos += melhorDelta;
            lista_tabu[melhorVizinho][corAnterior] = iteracao + 7;
            if (conflitos < melhorConflito)
            {
                melhorConflito = conflitos;
            }
        }
        iteracao++;
    }

    auto fim = high_resolution_clock::now();
    auto duracao = duration_cast<milliseconds>(fim - inicio).count();

    if (conflitos <= 0)
    {
        cout << "Solucao otima encontrada com " << n_cores << " cores" << endl;
    }
    else
    {
        cout << "Solucao nao otima encontrada" << endl;
    }

    cout << "Tempo de execucao: " << duracao << "ms ";
    cout << "/ Iterações: " << iteracao << endl;

    return conflitos <= 0;
};
#endif
