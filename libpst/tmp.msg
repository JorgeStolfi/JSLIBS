Seguem em anexo as figuras POV-Ray do gauge com baffle.

Para cortar o veludo voce precisa de um desenho do forro desenrolado.

Sejam
  dmin = diametro da base da forminha
  dmax = diametro do topo da forminha
  alt = altura da forminha
  
Calcule
  dr = (dmax - dmin)/2
  lado = sqrt(dr*dr + alt*alt)
  ang = 360*dr/lado
  rint = 0.5*dmin*lado/dr
  
O forro desenrolado é um pedaço de um anel,
com raio interno {rint}, raio externo {rint+lado},
e cobrindo um ângulo de {ang} graus visto do centro.

Vou tentar mandr um Postscript que faz as contas e plota.

