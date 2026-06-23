"""
Trabalho da 3ª Unidade – Grafos
Algoritmo: Fluxo Máximo (Edmonds-Karp) para distribuição de bolsas de sangue
Autores: Leonardo Machado Moreira e Sérgio Manuel Alves Moreira
UFRN – 2026
"""

from collections import deque

# =============================================================================
# ESTRUTURA DO GRAFO (lista de adjacência com capacidades)
# =============================================================================

class GrafoFluxo:
    """
    Grafo direcionado com capacidades para modelar a rede de fluxo.
    Representação: dicionário de dicionários  { u: { v: capacidade } }
    ex: aresta u → v com capacidade 5 é representada como grafo[u][v] = 5
    """

    def __init__(self):
        self.grafo = {}   # grafo residual
        self.nos   = set()

    def adicionar_aresta(self, u, v, capacidade):
        self.nos.add(u)
        self.nos.add(v)

        # aresta direta
        if u not in self.grafo:
            self.grafo[u] = {}
        self.grafo[u][v] = self.grafo[u].get(v, 0) + capacidade

        # aresta reversa (capacidade 0 para início)
        if v not in self.grafo:
            self.grafo[v] = {}
        if u not in self.grafo[v]:
            self.grafo[v][u] = 0

    def bfs(self, fonte, sumidouro, pai):
        """
        BFS no grafo residual – encontra caminho aumentante.
        Retorna True se existe caminho de fonte até sumidouro.
        """
        visitados = {fonte}
        fila = deque([fonte])

        while fila:
            u = fila.popleft()
            for v, capacidade in self.grafo[u].items():
                if v not in visitados and capacidade > 0:
                    visitados.add(v)
                    pai[v] = u # Se v for alcançado, u é seu pai no caminho
                    if v == sumidouro:
                        return True
                    fila.append(v)
        return False

    def edmonds_karp(self, fonte, sumidouro):
        """
        Algoritmo de Edmonds-Karp (Ford-Fulkerson com BFS).
        Retorna o valor do fluxo máximo.
        Complexidade: O(V * E²)
        """
        # guarda capacidades originais para depois calcular fluxo por aresta
        capacidade_original = {u: dict(vs) for u, vs in self.grafo.items()}

        fluxo_total = 0
        caminhos    = []   # guarda caminhos para exibição didática

        while True:
            pai = {}
            if not self.bfs(fonte, sumidouro, pai):
                break  # não há mais caminhos aumentantes

            # determina gargalo do caminho encontrado
            gargalo = float('inf')
            v = sumidouro
            while v != fonte:
                u = pai[v]
                gargalo = min(gargalo, self.grafo[u][v])
                v = u

            # reconstrói o caminho para exibição
            caminho = []
            v = sumidouro
            while v != fonte:
                caminho.append(v)
                v = pai[v]
            caminho.append(fonte)
            caminho.reverse()
            caminhos.append((caminho, gargalo))

            # atualiza capacidades residuais
            v = sumidouro
            while v != fonte:
                u = pai[v]
                self.grafo[u][v] -= gargalo
                self.grafo[v][u] += gargalo
                v = u

            fluxo_total += gargalo

        # fluxo realizado em cada aresta original = capacidade_original - capacidade_residual
        fluxo_por_aresta = {}
        for u, vizinhos in capacidade_original.items():
            for v, cap_orig in vizinhos.items():
                if cap_orig > 0:  # só arestas "reais", não as reversas artificiais
                    usado = cap_orig - self.grafo[u][v]
                    if usado > 0:
                        fluxo_por_aresta[(u, v)] = usado

        return fluxo_total, caminhos, fluxo_por_aresta

    def corte_minimo(self, fonte):
        """
        Após rodar edmonds_karp, identifica o corte mínimo:
        - Lado S: nós ainda alcançáveis pela fonte no grafo residual final
        - Lado T: os demais nós
        Retorna as arestas saturadas que cruzam de S para T (o gargalo).
        """
        alcancaveis = {fonte}
        fila = deque([fonte])
        while fila:
            u = fila.popleft()
            for v, cap in self.grafo[u].items():
                if cap > 0 and v not in alcancaveis:
                    alcancaveis.add(v)
                    fila.append(v)

        arestas_corte = []
        for u in alcancaveis:
            for v in self.grafo.get(u, {}):
                if v not in alcancaveis:
                    arestas_corte.append((u, v))

        return alcancaveis, arestas_corte


# =============================================================================
# DADOS DO PROBLEMA
# Estoque: baseado em proporções reais do perfil sanguíneo brasileiro
# (IBGE/Ministério da Saúde – tipo O+ é o mais comum: ~36% da população)
# =============================================================================

ESTOQUE_BOLSAS = {
    "O-":  20,
    "O+":  80,
    "A-":  15,
    "A+":  70,
    "B-":  10,
    "B+":  35,
    "AB-":  5,
    "AB+": 15,
}

DEMANDA_PACIENTES = {
    "Pac_O-":  12,
    "Pac_O+":  60,
    "Pac_A-":  10,
    "Pac_A+":  55,
    "Pac_B-":   8,
    "Pac_B+":  28,
    "Pac_AB-":  4,
    "Pac_AB+": 20,
}

# Cenário 2: escassez (ex.: período de baixa doação / alta demanda sazonal)
# Estoque de O- e O+ reduzido drasticamente, simulando falta do doador universal
ESTOQUE_ESCASSEZ = {
    "O-":   3,
    "O+":  15,
    "A-":  15,
    "A+":  70,
    "B-":  10,
    "B+":  35,
    "AB-":  5,
    "AB+": 15,
}

# Tabela de compatibilidade ABO/Rh
# chave = tipo do doador, valor = lista de receptores compatíveis
COMPATIBILIDADE = {
    "O-":  ["Pac_O-", "Pac_O+", "Pac_A-", "Pac_A+", "Pac_B-", "Pac_B+", "Pac_AB-", "Pac_AB+"],
    "O+":  ["Pac_O+", "Pac_A+", "Pac_B+", "Pac_AB+"],
    "A-":  ["Pac_A-", "Pac_A+", "Pac_AB-", "Pac_AB+"],
    "A+":  ["Pac_A+", "Pac_AB+"],
    "B-":  ["Pac_B-", "Pac_B+", "Pac_AB-", "Pac_AB+"],
    "B+":  ["Pac_B+", "Pac_AB+"],
    "AB-": ["Pac_AB-", "Pac_AB+"],
    "AB+": ["Pac_AB+"],
}

CAPACIDADE_INF = 10_000  # representa capacidade "infinita" nas arestas de compatibilidade


# =============================================================================
# CONSTRUÇÃO DO GRAFO DE FLUXO
# =============================================================================

def construir_grafo(estoque=None):
    """
    Monta a rede de fluxo em 4 camadas:
    Fonte (S) → Nós de oferta → Nós de demanda → Sumidouro (T)
    """
    if estoque is None:
        estoque = ESTOQUE_BOLSAS

    g = GrafoFluxo()
    FONTE     = "S"
    SUMIDOURO = "T"

    # Camada 1: Fonte → Oferta (capacidade = estoque disponível)
    for tipo, qtd in estoque.items():
        g.adicionar_aresta(FONTE, tipo, qtd)

    # Camada 2: Oferta → Demanda (capacidade infinita, limitada só pela compatibilidade)
    for doador, receptores in COMPATIBILIDADE.items():
        for receptor in receptores:
            g.adicionar_aresta(doador, receptor, CAPACIDADE_INF)

    # Camada 3: Demanda → Sumidouro (capacidade = necessidade do paciente)
    for tipo_pac, qtd in DEMANDA_PACIENTES.items():
        g.adicionar_aresta(tipo_pac, SUMIDOURO, qtd)

    return g, FONTE, SUMIDOURO


# =============================================================================
# EXIBIÇÃO DOS RESULTADOS
# =============================================================================

def exibir_resultados(estoque, fluxo_maximo, caminhos, fluxo_por_aresta, grafo, fonte, titulo):
    demanda_total = sum(DEMANDA_PACIENTES.values())

    print("=" * 64)
    print(f"  {titulo}")
    print("=" * 64)

    print("\n[ ESTOQUE DISPONÍVEL ]")
    for tipo, qtd in estoque.items():
        print(f"  {tipo:<5}: {qtd} bolsas")

    print("\n[ DEMANDA DOS PACIENTES ]")
    for tipo, qtd in DEMANDA_PACIENTES.items():
        print(f"  {tipo.replace('Pac_', ''):<5}: {qtd} pacientes")

    print(f"\n[ CAMINHOS AUMENTANTES ENCONTRADOS ({len(caminhos)}) ]")
    for i, (caminho, gargalo) in enumerate(caminhos, 1):
        print(f"  [{i:02d}] {' → '.join(caminho)}  |  fluxo: {gargalo}")

    print("\n[ RELATÓRIO DE ALOCAÇÃO (bolsas enviadas por rota) ]")
    alocacao_oferta_demanda = {
        (u, v): f for (u, v), f in fluxo_por_aresta.items()
        if u != fonte and v != "T"
    }
    for (u, v), f in sorted(alocacao_oferta_demanda.items()):
        print(f"  {u:<5} → {v.replace('Pac_', ''):<5} : {f} bolsas")

    print("\n[ RESULTADO FINAL ]")
    print(f"  Demanda total           : {demanda_total} pacientes")
    print(f"  Fluxo máximo            : {fluxo_maximo} pacientes atendidos")
    print(f"  Pacientes não atendidos : {demanda_total - fluxo_maximo}")

    if fluxo_maximo < demanda_total:
        print("\n  ⚠ Estoque insuficiente para atender toda a demanda.")
        alcancaveis, arestas_corte = grafo.corte_minimo(fonte)
        print("  → Corte mínimo (gargalo do sistema):")
        for u, v in arestas_corte:
            if u == fonte:
                print(f"     Estoque insuficiente de: {v}")
    else:
        print("\n  ✓ Toda a demanda pode ser atendida com o estoque atual.")

    print("=" * 64)


# =============================================================================
# VISUALIZAÇÃO DO GRAFO
# =============================================================================

def desenhar_grafo(fluxo_por_aresta, estoque, arquivo_saida):
    """
    Desenha a rede de fluxo em 4 camadas usando matplotlib,
    destacando as arestas que carregam fluxo > 0.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(13, 8))

    tipos = list(estoque.keys())
    pacientes = list(DEMANDA_PACIENTES.keys())

    pos = {"S": (0, len(tipos) / 2)}
    for i, t in enumerate(tipos):
        pos[t] = (1, len(tipos) - i)
    for i, p in enumerate(pacientes):
        pos[p] = (2, len(pacientes) - i)
    pos["T"] = (3, len(pacientes) / 2)

    # desenha todas as arestas de compatibilidade em cinza claro (fundo)
    for doador, receptores in COMPATIBILIDADE.items():
        for receptor in receptores:
            x1, y1 = pos[doador]
            x2, y2 = pos[receptor]
            ax.plot([x1, x2], [y1, y2], color="#dddddd", linewidth=0.6, zorder=1)

    # destaca arestas com fluxo realizado
    for (u, v), f in fluxo_por_aresta.items():
        if u in pos and v in pos and u != "S" and v != "T":
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            ax.plot([x1, x2], [y1, y2], color="#c0392b", linewidth=1.2 + f / 15, zorder=2, alpha=0.85)

    # nós S e T
    for n in ["S", "T"]:
        x, y = pos[n]
        ax.scatter(x, y, s=900, color="#2e4d7b", zorder=3)
        ax.text(x, y, n, ha="center", va="center", color="white", fontsize=11, fontweight="bold", zorder=4)

    # nós de oferta
    for t in tipos:
        x, y = pos[t]
        ax.scatter(x, y, s=1300, color="#a9cce3", edgecolor="#2e4d7b", zorder=3)
        ax.text(x, y, f"{t}\n{estoque[t]}", ha="center", va="center", fontsize=8, zorder=4)

    # nós de demanda
    for p in pacientes:
        x, y = pos[p]
        label = p.replace("Pac_", "")
        ax.scatter(x, y, s=1300, color="#f5b7b1", edgecolor="#922b21", zorder=3)
        ax.text(x, y, f"{label}\n{DEMANDA_PACIENTES[p]}", ha="center", va="center", fontsize=8, zorder=4)

    legenda = [
        Line2D([0], [0], color="#c0392b", lw=2, label="Fluxo realizado (bolsas alocadas)"),
        Line2D([0], [0], color="#dddddd", lw=2, label="Compatibilidade possível (sem fluxo)"),
        mpatches.Patch(color="#a9cce3", label="Tipo sanguíneo (oferta)"),
        mpatches.Patch(color="#f5b7b1", label="Grupo de pacientes (demanda)"),
    ]
    ax.legend(handles=legenda, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=2, fontsize=9)

    ax.set_title("Rede de Fluxo Máximo — Distribuição de Bolsas de Sangue", fontsize=13, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(arquivo_saida, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n[Gráfico salvo em: {arquivo_saida}]")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # ---------- Cenário 1: estoque normal ----------
    grafo1, fonte, sumidouro = construir_grafo(ESTOQUE_BOLSAS)
    fluxo1, caminhos1, alocacao1 = grafo1.edmonds_karp(fonte, sumidouro)
    exibir_resultados(ESTOQUE_BOLSAS, fluxo1, caminhos1, alocacao1, grafo1, fonte,
                       "CENÁRIO 1 — ESTOQUE NORMAL")
    desenhar_grafo(alocacao1, ESTOQUE_BOLSAS, "grafo_cenario1_normal.png")

    print("\n\n")

    # ---------- Cenário 2: escassez ----------
    grafo2, fonte, sumidouro = construir_grafo(ESTOQUE_ESCASSEZ)
    fluxo2, caminhos2, alocacao2 = grafo2.edmonds_karp(fonte, sumidouro)
    exibir_resultados(ESTOQUE_ESCASSEZ, fluxo2, caminhos2, alocacao2, grafo2, fonte,
                       "CENÁRIO 2 — ESCASSEZ DE TIPO O (doador universal)")
    desenhar_grafo(alocacao2, ESTOQUE_ESCASSEZ, "grafo_cenario2_escassez.png")
