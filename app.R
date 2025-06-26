library(GGally)
library(dplyr)



indice_inversion_vecteurs <- function(v1, v2) {
  if (length(v1) != length(v2)) {
    stop("Les deux vecteurs doivent avoir la même longueur.")
  }
  
  r1 <- rank(-v1, ties.method = "min")
  r2 <- rank(-v2, ties.method = "min")
  
  indice <- -cor(r1, r2, method = "spearman")
  return(indice)
}


# On extrait les données DOWN et UP dans deux vecteurs nommés
df_down <- score_XR_pur %>%
  filter(BaseName_pur == "XBP1_pur_DOWN") %>%
  select(Sample_pur, Score_pur) %>%
  rename(DOWN = Score_pur)

df_up <- score_XR_pur %>%
  filter(BaseName_pur == "XBP1_pur_UP") %>%
  select(Sample_pur, Score_pur) %>%
  rename(UP = Score_pur)

# Fusion des deux
df_scores <- left_join(df_down, df_up, by = "Sample_pur")

# Rang (1 = plus haut score)
df_scores <- df_scores %>%
  mutate(
    rang_DOWN = rank(-DOWN, ties.method = "min"),
    rang_UP = rank(-UP, ties.method = "min")
  )

# Calcul des seuils de déciles
seuils_down <- quantile(df_scores$DOWN, probs = seq(0.1, 1, 0.1), na.rm = TRUE)
seuils_up <- quantile(df_scores$UP, probs = seq(0.1, 1, 0.1), na.rm = TRUE)

# Attribution du décile à chaque individu
df_scores <- df_scores %>%
  mutate(
    decile_DOWN = cut(DOWN, breaks = c(-Inf, seuils_down), labels = 1:10, include.lowest = TRUE),
    decile_UP = cut(UP, breaks = c(-Inf, seuils_up), labels = 1:10, include.lowest = TRUE)
  )

# Trouver les individus les plus proches de chaque seuil de décile
get_individus_deciles <- function(scores, seuils) {
  sapply(seuils, function(seuil) {
    which.min(abs(scores - seuil))
  })
}

ind_down <- get_individus_deciles(df_scores$DOWN, seuils_down)
ind_up <- get_individus_deciles(df_scores$UP, seuils_up)

individus_deciles_down <- df_scores$Sample_pur[ind_down]
individus_deciles_up <- df_scores$Sample_pur[ind_up]

# Afficher les résultats
print("Rangs et déciles par échantillon :")
print(df_scores)

print("Individus proches des seuils de déciles (DOWN) :")
print(data.frame(decile = 1:10, Sample_pur = individus_deciles_down, seuil = round(seuils_down, 3)))

print("Individus proches des seuils de déciles (UP) :")
print(data.frame(decile = 1:10, Sample_pur = individus_deciles_up, seuil = round(seuils_up, 3)))


# Récupérer les rangs des individus déciles (pas mêmes individus dans DOWN et UP)
individus_deciles <- df_scores %>%
  filter(Sample_pur %in% individus_deciles_up) %>%  # ou individus_deciles_up, c’est équivalent
  select(Sample_pur, rang_DOWN, rang_UP)

# Calcul de l’indice d’inversion (corrélation de Spearman inversée)
indice_inversion_deciles <- -cor(individus_deciles$rang_DOWN, individus_deciles$rang_UP, method = "spearman")

print(paste("Indice d'inversion des rangs pour les déciles :", round(indice_inversion_deciles, 3)))


# Préparer le data.frame pour ggparcoord
df_parcoord <- individus_deciles %>%
  mutate(Sample = Sample_pur) %>%
  select(Sample, rang_DOWN, rang_UP)

# Tracé parallèle
ggparcoord(df_parcoord,
           columns = 2:3,  # rang_DOWN et rang_UP
           groupColumn = 1,  # Sample
           showPoints = TRUE,
           scale = "globalminmax") +
  
  theme_minimal() +
  labs(
    title = paste0("Inversion des rangs des déciles pour les déciles UP— Indice : ", round(indice_inversion_deciles, 2)),
    x = "Condition (DOWN → UP)",
    y = "Rang (1 = plus exprimé)"
  )



scores_deciles_DOWN <- df_scores %>%
  filter(Sample_pur %in% individus_deciles_down) %>%
  pull(DOWN)

scores_deciles_UP <- df_scores %>%
  filter(Sample_pur %in% individus_deciles_down) %>%
  pull(UP)

# Calcul de l’indice d’inversion des rangs des déciles :
indice_d <- indice_inversion_vecteurs(scores_deciles_DOWN, scores_deciles_UP)
print(paste("Indice d'inversion des rangs des déciles :", round(indice_d, 3)))

