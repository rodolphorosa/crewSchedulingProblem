//============================================================================
// Name        : TopicosEspeciaisOtimizacao.cpp
// Author      : Rodolpho Rosa da Silva
// Version     : 1.0
// Copyright   : 
// Description : Iterated Local Search with Variable Neighbourhood Descent as local search method to solve the Crew Scheduling Problem
//============================================================================

#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

const int MAX_TEMPO = 480;
const int NIVEL_MIN = 2;
const int NIVEL_MAX = 10;
const int MAX_ITERACOES = 1000;

class ArcoTransicao {
	int origem;
	int destino;
	int custo;
	public:
		ArcoTransicao(int origem, int destino, int custo): origem(origem), destino(destino), custo(custo) {}

		int retornarOrigem() { return origem; }

		int retornarDestino() { return destino; }

		int retornarCusto() { return custo; }

		bool temDestino(int _destino) {	return destino == _destino; }
};

class Tarefa {
	int tarefa;
	int inicio;
	int termino;
	vector<ArcoTransicao> transicoes;
	public:
		Tarefa(int tarefa, int inicio, int termino): tarefa(tarefa), inicio(inicio), termino(termino) {}

		int retornarTempo() { return termino - inicio; }
		int retornarInicio() { return inicio; }
		int retornarTermino() { return termino; }
		void adicionaArcoTransicao(ArcoTransicao arco) { transicoes.push_back(arco); }
		int retornarNome() { return tarefa; }
		int retornarNumeroTransicoes() { return (int) transicoes.size(); }

		int retornarCustoTransicao(int tarefa) {
			for(int i=0; i<retornarNumeroTransicoes(); i++) {
				if(transicoes[i].temDestino(tarefa)) return transicoes[i].retornarCusto();
			}
			return 0;
		}

		bool temTransicao(int tarefa) {
			if(tarefa == 0) return true;
			for(int i=0; i<retornarNumeroTransicoes(); i++) {
				if(transicoes[i].temDestino(tarefa)) return true;
			}
			return false;
		}

		bool equals(Tarefa t) {
			return this->tarefa == t.retornarNome();
		}

		vector<ArcoTransicao> retornarTransicoes() { return transicoes; }
};

class Jornada {
	vector<Tarefa> tarefas;
	int custo;
	int tempo;
	int max_tempo;
	public:
		Jornada(int max_tempo): max_tempo(max_tempo) { custo = 0; tempo = 0; }

		void adicionarTarefa(Tarefa t) {
			if(tarefas.empty()) {
				tarefas.push_back(t);
				tempo += t.retornarTempo();
			} else {
				Tarefa ultima = tarefas.back();
				tarefas.push_back(t);
				custo = calcularCusto();
				tempo = calcularTempo();
			}
		}

		bool jornadaViavel(Tarefa t, int max_tempo) {
			if(tarefas.empty()) return true;
			Tarefa ultima = tarefas.back();
			return (ultima.temTransicao(t.retornarNome()))
					&& (tempo + t.retornarTempo() + (t.retornarInicio() - ultima.retornarTermino())) <= max_tempo
					&& ultima.retornarTermino() <= t.retornarInicio();
		}

		void realizarTroca(Tarefa nova, Tarefa atual) {
			for(int i=0; i<(int)tarefas.size(); i++) {
				if(tarefas[i].retornarNome() == atual.retornarNome()) {
					tarefas[i] = nova;
					custo = calcularCusto();
					tempo = calcularTempo();
				}
			}
		}

		void alocarTarefa(Tarefa nova, int posicao) {
			vector<Tarefa>::iterator iterador = tarefas.begin()+posicao;
			tarefas.insert(iterador, nova);
			custo = calcularCusto();
			tempo = calcularTempo();
		}

		void removerTarefa(Tarefa tarefa) {
			int pos;
			for(pos=0; pos<(int)tarefas.size(); pos++) {
				if(tarefas[pos].retornarNome() == tarefa.retornarNome()) {
					tarefas.erase(tarefas.begin()+pos);
					tempo = calcularTempo();
					custo = calcularCusto();
				}
			}
		}

		bool restricoesVioladas() {
			if(tempo > max_tempo) return true;
			for(int i=0; i<(int)(tarefas.size())-1; i++) {
				if(!tarefas[i].temTransicao(tarefas[i+1].retornarNome())) return true;
			}
			return false;
		}

		int retornarTempo() { return tempo; }
		int retornarCusto() { return custo; }
		vector<Tarefa> retornarTarefas() { return tarefas; }

	private:
		int calcularCusto() {
			int soma = 0;
			for(int i=0; i<(int)tarefas.size()-1; i++) {
				soma += tarefas[i].retornarCustoTransicao(tarefas[i+1].retornarNome());
			}
			return soma;
		}

		int calcularTempo() {
			int soma = tarefas[0].retornarTempo();
			for(int i=1; i<(int)tarefas.size(); i++) {
				soma += tarefas[i].retornarTempo() + (tarefas[i].retornarInicio() - tarefas[i-1].retornarTermino());
			}
			return soma;
		}
};

vector<int> sortearJornadas(int max_tarefas);
int custoSolucao(vector<Jornada>jornadas);
void imprimirSolucao(vector<Jornada>jornadas);

vector<Jornada> gerarSolucaoInicial(vector<Tarefa>tarefas, int max_tempo);
vector<Jornada> buscarSolucaoVizinha(vector<Jornada>jornadas);
vector<Jornada> aplicarPerturbacao(vector<Jornada>solucaoBase, int nivel);
vector<Jornada> variableNeighhoodDescent(vector<Jornada>solucaoInicial);
vector<Jornada> aplicarCriterioAceitacao(vector<Jornada>solucaoCandidata, vector<Jornada>solucaoOtima, int &nivelPerturbacao, int nivelMIn, int nivelMax);
vector<Jornada> iteratedLocalSearch(vector<Jornada>solucaoInicial, int nivelMin, int nivelMax);

int main() {

	srand(time(NULL));

	ifstream arquivo;
	arquivo.open("./dataset/csp50.txt");

	int qtde_tarefas;
	int max_tempo;

	arquivo >> qtde_tarefas >> max_tempo;

	vector<Tarefa> tarefas;
	vector<Jornada> jornadas;

	for(int i=0; i<qtde_tarefas; i++) {
		int inicio, termino;
		arquivo >> inicio >> termino;
		Tarefa t (i+1, inicio, termino);
		tarefas.push_back(t);
	}

	int origem, destino, custo;

	while(arquivo >> origem >> destino >> custo) {
		ArcoTransicao transicao (origem, destino, custo);
		tarefas[origem-1].adicionaArcoTransicao(transicao);
	}

	vector<Jornada>solucaoInicial = gerarSolucaoInicial(tarefas, max_tempo);

	//imprimirSolucao(solucaoInicial);

	cout << "Solucao inicial:\nnumero de jornadas: " << solucaoInicial.size() << "\ncusto da solucao inicial: " << custoSolucao(solucaoInicial) << "\n\n";

	vector<Jornada>solucaoOtima = iteratedLocalSearch(solucaoInicial, NIVEL_MIN, NIVEL_MAX);

	imprimirSolucao(solucaoOtima);

	cout << "Melhor solucao:\nnumero de jornadas: " << solucaoOtima.size() << "\nCusto da melhor solucao: " << custoSolucao(solucaoOtima) << "\n\n";

	cout << "Taxa de melhora: " << (float) ((custoSolucao(solucaoInicial) - custoSolucao(solucaoOtima)) * 100) / (custoSolucao(solucaoInicial)) << "\n";

	return 0;
}

vector<Jornada> gerarSolucaoInicial(vector<Tarefa>tarefas, int tempo_maximo) {
	int numero_tarefas = (int) tarefas.size();
	vector<Jornada>jornadas;

	Jornada j(tempo_maximo);
	j.adicionarTarefa(tarefas[0]);
	jornadas.push_back(j);

	for(int i=1; i<numero_tarefas; i++) {
		int jornada = -1;

		for(int j=0; j<(int)jornadas.size(); j++) {
			if(jornadas[j].jornadaViavel(tarefas[i], tempo_maximo)) {
				jornada = j;
			}
		}

		if(jornada == -1) {
			Jornada j(tempo_maximo);
			j.adicionarTarefa(tarefas[i]);
			jornadas.push_back(j);
		}
		else { jornadas[jornada].adicionarTarefa(tarefas[i]); }
	}
	return jornadas;
}

vector<Jornada> buscarSolucaoVizinha(vector<Jornada>jornadas) {

	while(true) {
		vector<int>n = sortearJornadas((int)jornadas.size());
		vector<int>ts;

		/* Impossibilita troca entre jornadas iguais */
		if(n[0] == n[1]) continue;

		Jornada j = jornadas[n[0]];
		Jornada j_linha = jornadas[n[1]];

		for(int i=0; i<2; i++) {
			int random = (int)rand() % jornadas[n[i]].retornarTarefas().size();
			ts.push_back(random);
		}

		Tarefa t = j.retornarTarefas()[ts[0]];
		Tarefa t_linha = j_linha.retornarTarefas()[ts[1]];

		/* Verifica se ha uma realocacao valida */
		for(int i=0; i<=(int)j_linha.retornarTarefas().size(); i++) {
			Jornada aux = j;
			Jornada aux_linha = j_linha;
			aux.removerTarefa(t);
			aux_linha.alocarTarefa(t, i);
			if(!aux.restricoesVioladas() && !aux_linha.restricoesVioladas()) {
				if(aux.retornarTarefas().size() == 0) {
					jornadas[n[1]] = aux_linha;
					jornadas.erase(jornadas.begin()+n[0]);
				} else {
					jornadas[n[0]] = aux;
					jornadas[n[1]] = aux_linha;
				}
				return jornadas;
			}
		}

		/* Descarta trocas simetricas */
		if((int)j.retornarTarefas().size() == 1 && (int)j_linha.retornarTarefas().size() == 1) continue;

		j = jornadas[n[0]];
		j_linha = jornadas[n[1]];

		/* Realiza a troca entre duas tarefas */
		j.realizarTroca(t_linha, t);
		j_linha.realizarTroca(t, t_linha);

		if(!j.restricoesVioladas() && !j_linha.restricoesVioladas()) {
			jornadas[n[0]] = j;
			jornadas[n[1]] = j_linha;
			return jornadas;
		}
	}
}

vector<Jornada> aplicarPerturbacao(vector<Jornada>solucaoBase, int nivel) {

	vector<Jornada>solucao = solucaoBase;
	vector<Jornada>solucaoOtima;

	for(int i=0; i<nivel; i++) {
		solucaoOtima = buscarSolucaoVizinha(solucao);
		solucao = solucaoOtima;
	}
	return solucao;
}

vector<Jornada>variableNeighhoodDescent(vector<Jornada>solucaoInicial) {
	vector<Jornada>solucaoOtima = solucaoInicial;
	int r = 2;
	int k = 0;
	while(k < r) {
		vector<Jornada> solucao = aplicarPerturbacao(solucaoOtima, k);
		if(custoSolucao(solucao) < custoSolucao(solucaoOtima)) {
			solucaoOtima = solucao;
			k = 0;
		} else {
			k = k + 1;
		}
	}
	return solucaoOtima;
}

vector<Jornada> aplicarCriterioAceitacao(vector<Jornada>solucaoCandidata, vector<Jornada>solucaoOtima, int &nivelPerturbacao, int nivelMin, int nivelMax) {
	vector<Jornada>solucao;

	if(custoSolucao(solucaoCandidata) < custoSolucao(solucaoOtima)) {
		nivelPerturbacao = nivelMin;
		solucao = solucaoCandidata;
	} else {
		nivelPerturbacao = std::min(nivelPerturbacao+1, nivelMax);
		solucao = solucaoOtima;
	}
	return solucao;
}

vector<Jornada> iteratedLocalSearch(vector<Jornada>solucaoInicial, int nivelMin, int nivelMax) {
	vector<Jornada>solucaoOtima = solucaoInicial;
	int nivelPerturbacao;
	int iteracoes = 0;

	solucaoOtima = variableNeighhoodDescent(solucaoInicial);
	nivelPerturbacao = nivelMin;

	while(true) {
		vector<Jornada>solucaoLinha = aplicarPerturbacao(solucaoOtima, nivelPerturbacao);
		vector<Jornada>solucaoOtimaLinha = variableNeighhoodDescent(solucaoLinha);
		solucaoOtima = aplicarCriterioAceitacao(solucaoOtima, solucaoOtimaLinha, nivelPerturbacao, nivelMin, nivelMax);
		if(iteracoes > MAX_ITERACOES) return solucaoOtima;
		iteracoes++;
	}
}

vector<int> sortearJornadas(int max_jornadas) {
	vector<int>jornadas;
	for(int i=0; i<2; i++) {
		int n = (int) rand() % max_jornadas;
		jornadas.push_back(n);
	}
	return jornadas;
}

void imprimirSolucao(vector<Jornada>jornadas) {
	for(int i=0; i<(int)jornadas.size(); i++) {
		Jornada j = jornadas[i];
		cout << i << " : ";
		for(int k=0; k<(int)j.retornarTarefas().size(); k++) {
			cout << j.retornarTarefas()[k].retornarNome() << " ";
		}
		cout << " custo : " << j.retornarCusto() << " tempo: " << j.retornarTempo() << "\n";
	}
}

int custoSolucao(vector<Jornada>jornadas) {
	int soma = 0;
	for(int i=0; i<(int)jornadas.size(); i++) {
		soma = soma + jornadas[i].retornarCusto();
	}
	return soma;
}
