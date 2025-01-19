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

struct ResponseMaxColor
{
    vector<int> cores;
    int maxColor;
};

struct ResponseBuscaLocal
{
    vector<int> solucao_encontrada;
    int num_conflitos;
};

class Problem
{
public:
    vector<int> calcClique(Graph &graph);
    vector<int> troca_cores_aleatorio(vector<int> &cores, Graph &graph);
    vector<int> troca(vector<int> &cores, int pos1, int pos2);
    vector<int> substitui_cor_menos_frequente(vector<int> &cores, Graph &graph);
    ResponseMaxColor getMaxColor(Graph &Graph, int num_vertices, int &grauMaximo);
    HeuristicaConstrutivaResponse heuristicaConstrutiva(Graph &graph, int k);
    int calcConflitoVertice(vector<int> &cores, Graph &graph, int vertice);
    int calcConflitoTotal(vector<int> &cores, Graph &graph);
    vector<int> calcConflitoGraph(vector<int> &cores, Graph &graph);
    vector<int> ordenaVerticesPorConflito(vector<int> &solucao, Graph &graph);
    unordered_map<int, int> calcFrequenciaCores(vector<int> &solucao);
    vector<int> rankCorPorMenorFrequencia(unordered_map<int, int> &frequencias);
    ResponseBuscaLocal busca_local_melhor_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes);
    ResponseBuscaLocal busca_local_primeira_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes);
    ResponseBuscaLocal busca_local_primeira_melhora_cor_menos_frequente(Graph &graph, vector<int> &solucao_inicial, int iteracoes);
};

vector<int> Problem::rankCorPorMenorFrequencia(unordered_map<int, int> &frequencias)
{
    unordered_map<int, int> frequencias_aux = frequencias;
    vector<int> rank_cores;

    while (!frequencias_aux.empty())
    {
        auto menorFreqIt = min_element(
            frequencias.begin(), frequencias.end(),
            [](const pair<int, int> &a, const pair<int, int> &b)
            {
                return a.second < b.second;
            });

        rank_cores.push_back(menorFreqIt->first);

        frequencias.erase(menorFreqIt);
    }
    return rank_cores;
}

ResponseBuscaLocal Problem::busca_local_primeira_melhora_cor_menos_frequente(Graph &graph, vector<int> &solucao_inicial, int iteracoes)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;
    int melhor_custo = *max_element(solucao.begin(), solucao.end());
    int melhor_conflito = calcConflitoTotal(solucao, graph);
    vector<int> vertices_ordenados_conflitos = ordenaVerticesPorConflito(solucao, graph);
    unordered_map<int, int> frequenciaCores = calcFrequenciaCores(melhor_vizinho);
    vector<int> rank_cores_menor_frequencia = rankCorPorMenorFrequencia(frequenciaCores);

    int cont_iteracoes = 0;
    for (int i = 0; i < vertices_ordenados_conflitos.size(); i++)
    {

        for (CoupleVertice aresta : graph.get_neighbors(i))
        {

            int vertice1 = i;
            int vertice2 = aresta.vertice2;
            vector<int> vizinho = troca(melhor_vizinho, vertice1, vertice2);
            int conflito_vizinho = calcConflitoTotal(vizinho, graph);
            if (conflito_vizinho < melhor_conflito)
            {

                melhor_vizinho = vizinho;
                melhor_conflito = conflito_vizinho;
                break;
            }
            cont_iteracoes++;
            if (cont_iteracoes == iteracoes)
            {
                break;
            }
        }
        if (cont_iteracoes == iteracoes)
        {
            break;
        }
    }

    ResponseBuscaLocal response;
    response.solucao_encontrada = melhor_vizinho;
    response.num_conflitos = melhor_conflito;

    return response;
}

vector<int> Problem::substitui_cor_menos_frequente(vector<int> &cores, Graph &graph)
{
    int num_vertices = graph.get_num_vertices();
    unordered_map<int, int> frequencias = calcFrequenciaCores(cores);
    int cor_menos_frequente = -1;
    int menor_frequencia = num_vertices + 1;
    for (auto &par : frequencias)
    {
        if (par.second < menor_frequencia)
        {
            menor_frequencia = par.second;
            cor_menos_frequente = par.first;
        }
    }

    for (int i = 0; i < num_vertices; ++i)
    {
        if (cores[i] == cor_menos_frequente)
        {
            set<int> cores_vizinhas;
            for (CoupleVertice vizinho : graph.get_neighbors(i))
            {
                cores_vizinhas.insert(cores[vizinho.vertice2]);
            }

            for (int nova_cor = 1;; ++nova_cor)
            {
                if (cores_vizinhas.find(nova_cor) == cores_vizinhas.end())
                {
                    cores[i] = nova_cor;
                    break;
                }
            }
        }
    }

    return cores;
};

unordered_map<int, int> Problem::calcFrequenciaCores(vector<int> &solucao)
{
    unordered_map<int, int> frequencias;

    for (int cor : solucao)
    {
        frequencias[cor]++;
    }

    return frequencias;
};

vector<int> Problem::ordenaVerticesPorConflito(vector<int> &solucao, Graph &graph)
{

    vector<int> conflitos = calcConflitoGraph(solucao, graph);

    vector<int> vertices_ordenados_conflitos;

    int maior = INT32_MAX;
    int maior_atual = 0;
    int maior_vertice_atual = 0;
    int pos_ordenacao = 0;

    for (int j = 0; j < graph.get_num_vertices(); j++)
    {
        for (int i = 0; i < conflitos.size(); i++)
        {
            if (conflitos[i] > maior_atual && conflitos[i] < maior)
            {
                maior_atual = conflitos[i];
                maior_vertice_atual = i;
            }
        }

        maior = maior_atual;
        vertices_ordenados_conflitos.push_back(maior_vertice_atual);
    }

    return vertices_ordenados_conflitos;
}

ResponseBuscaLocal Problem::busca_local_primeira_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;
    // int melhor_custo = *max_element(solucao.begin(), solucao.end());
    int melhor_conflito = calcConflitoTotal(solucao, graph);
    vector<int> vertices_ordenados_conflitos = ordenaVerticesPorConflito(solucao, graph);

    int cont_iteracoes = 0;
    for (int i = 0; i < vertices_ordenados_conflitos.size(); i++)
    {

        for (CoupleVertice aresta : graph.get_neighbors(i))
        {

            int vertice1 = i;
            int vertice2 = aresta.vertice2;
            vector<int> vizinho = troca(melhor_vizinho, vertice1, vertice2);
            int conflito_vizinho = calcConflitoTotal(vizinho, graph);
            if (conflito_vizinho < melhor_conflito)
            {
                // cout << "feitas " << cont_iteracoes << " iteracoes  n conflitos " << conflito_vizinho << endl;
                melhor_vizinho = vizinho;
                melhor_conflito = conflito_vizinho;
                break;
            }
            cont_iteracoes++;
            if (cont_iteracoes == iteracoes)
            {
                break;
            }
        }
        if (cont_iteracoes == iteracoes)
        {
            break;
        }

        // vector<int> vizinho = troca_cores_aleatorio(melhor_vizinho, graph);

        // // int custo_vizinho = *max_element(vizinho.begin(), vizinho.end());
        // int conflito_vizinho = calcConflitoTotal(vizinho, graph);

        // if (conflito_vizinho < melhor_conflito)
        // {
        //     cout << "feitas " << i << " iteracoes " << endl;
        //     melhor_vizinho = vizinho;
        //     melhor_conflito = conflito_vizinho;
        // }
        // if (i == iteracoes - 1)
        // {
        //     cout << "feitas " << iteracoes - 1 << " iteracoes " << endl;
        // }
    }

    ResponseBuscaLocal response;
    response.solucao_encontrada = melhor_vizinho;
    response.num_conflitos = melhor_conflito;

    return response;
}

ResponseBuscaLocal Problem::busca_local_melhor_melhora(Graph &graph, vector<int> &solucao_inicial, int iteracoes)
{
    vector<int> solucao = solucao_inicial;
    vector<int> melhor_vizinho = solucao;
    int melhor_custo = *max_element(solucao.begin(), solucao.end());
    int melhor_conflito = calcConflitoTotal(solucao, graph);
    vector<int> vertices_ordenados_conflitos = ordenaVerticesPorConflito(solucao, graph);

    int cont_iteracoes = 0;

    for (int i = 0; i < vertices_ordenados_conflitos.size(); i++)
    {
        vector<int> melhor_vizinho_atual = melhor_vizinho;
        for (CoupleVertice aresta : graph.get_neighbors(i))
        {

            int vertice1 = i;
            int vertice2 = aresta.vertice2;
            vector<int> vizinho = troca(melhor_vizinho_atual, vertice1, vertice2);
            int conflito_vizinho = calcConflitoTotal(vizinho, graph);
            if (conflito_vizinho < melhor_conflito)
            {
                // cout << "feitas " << cont_iteracoes << " iteracoes  n conflitos " << conflito_vizinho << endl;
                melhor_vizinho = vizinho;
                melhor_conflito = conflito_vizinho;
            }
            cont_iteracoes++;
            if (cont_iteracoes == iteracoes)
            {
                break;
            }
        }

        // int custo_vizinho = *max_element(vizinho.begin(), vizinho.end());

        // if (i == iteracoes - 1)
        // {
        //     cout << "feitas " << iteracoes - 1 << " iteracoes " << endl;
        // }
        if (cont_iteracoes == iteracoes)
        {
            break;
        }
    }

    ResponseBuscaLocal response;
    response.solucao_encontrada = melhor_vizinho;
    response.num_conflitos = melhor_conflito;

    return response;
}

int Problem::calcConflitoTotal(vector<int> &cores, Graph &graph)
{
    int conflitos = 0;

    for (int vertice = 0; vertice < graph.get_num_vertices(); vertice++)
    {
        conflitos += calcConflitoVertice(cores, graph, vertice);
    }
    return conflitos;
}

vector<int> Problem::calcConflitoGraph(vector<int> &cores, Graph &graph)
{
    vector<int> conflitos(graph.get_num_vertices(), 0);

    for (int vertice = 0; vertice < graph.get_num_vertices(); vertice++)
    {
        conflitos[vertice] = calcConflitoVertice(cores, graph, vertice);
    }
    return conflitos;
}

int Problem::calcConflitoVertice(vector<int> &cores, Graph &graph, int vertice)
{
    int conflito = 0;
    for (CoupleVertice vizinho : graph.get_neighbors(vertice))
    {
        if (cores[vizinho.vertice2] == cores[vertice])
            conflito++;
    }
    return conflito;
}

vector<int> Problem::troca_cores_aleatorio(vector<int> &cores, Graph &graph)
{
    vector<int> cores_aux = cores;
    srand(time(0));
    int vertice1 = rand() % graph.get_num_vertices();
    int vertice2 = 0;
    do
    {
        vertice2 = rand() % graph.get_num_vertices();
    } while (vertice1 == vertice2);

    swap(cores_aux[vertice1], cores_aux[vertice2]);

    return cores_aux;
}

vector<int> Problem::troca(vector<int> &cores, int pos1, int pos2)
{
    vector<int> cores_aux = cores;
    swap(cores_aux[pos1], cores_aux[pos2]);
    return cores_aux;
}

vector<int> Problem::calcClique(Graph &graph)
{
    int n = graph.get_num_vertices();
    vector<int> vertices(n);
    for (int i = 0; i < n; ++i)
        vertices[i] = i;

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

ResponseMaxColor Problem::getMaxColor(Graph &graph, int num_vertices, int &grauMaximo)
{
    vector<int> cores(num_vertices, -1);
    vector<int> saturacao(num_vertices, 0);
    vector<int> graus(num_vertices, 0);

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
                if (saturacao[vertice] > maiorSaturacao ||
                    (saturacao[vertice] == maiorSaturacao && graus[vertice] > maiorGrau))
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

    ResponseMaxColor response;
    response.cores = cores;
    response.maxColor = maxCor;

    return response;
}

HeuristicaConstrutivaResponse Problem::heuristicaConstrutiva(Graph &graph, int initialK = 10)
{

    vector<int> clique = calcClique(graph);
    int k = clique.size();

    int num_arestas = 0;
    for (int i = 0; i < graph.get_num_vertices(); i++)
    {
        num_arestas += graph.get_neighbors(i).size();
    }

    int num_vertices = graph.get_num_vertices();
    int grauMaximo = 0;
    ResponseMaxColor responseMaxColor = getMaxColor(graph, num_vertices, grauMaximo);
    int maxCor = responseMaxColor.maxColor;

    int conflitos = calcConflitoTotal(responseMaxColor.cores, graph);

    cout << "conflito antes: " << conflitos << " maxCorAntes: " << maxCor << endl;

    // ResponseBuscaLocal responseBuscaLocal;
    // responseBuscaLocal.solucao_encontrada = responseMaxColor.cores;
    // responseBuscaLocal.num_conflitos = conflitos;

    // ResponseBuscaLocal busca = busca_local_primeira_melhora(graph, responseMaxColor.cores, 1000);
    // cout << "busca local primeira melhora -- conflito: " << busca.num_conflitos << endl;

    // busca = busca_local_melhor_melhora(graph, responseMaxColor.cores, 1000);
    // cout << "busca local melhor melhora -- conflito depois: " << busca.num_conflitos << endl;

    // for (int i = 0; i < 100; i++)
    // {
    //

    //     if (busca.num_conflitos > responseBuscaLocal.num_conflitos)
    //     {
    //         break;
    //     }

    //     responseBuscaLocal = busca;
    // }

    // int maxCorDepois = *max_element(responseBuscaLocal.solucao_encontrada.begin(), responseBuscaLocal.solucao_encontrada.end());
    // cout << "conflito depois: " << busca.num_conflitos << " maxCorDepois: " << maxCorDepois << endl;

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
#endif