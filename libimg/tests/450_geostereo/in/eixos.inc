// ======================================================================
// EIXOS DE COORDENADAS

#macro eixo(ponta,cor,nome)
  #local raio = 0.0075*vlength(ponta);
  union{
    sphere{ <0,0,0>, raio }
    cylinder{ <0,0,0>, 0.91*ponta, raio }
    cone{ 0.90*ponta, 3*raio, ponta, 0.00 }
    texture{ pigment{ color rgb cor } finish{ ambient 0.5 diffuse 0.5 } }
  }
#end
#macro eixos(tamanho)
  union{ 
    object{ eixo( <tamanho,0,0>, <1.0,0.2,0.2>, "X") }
    object{ eixo( <0,tamanho,0>, <0.0,0.8,0.0>, "Y") }
    object{ eixo( <0,0,tamanho>, <0.3,0.3,1.0>, "Z") }
  }
#end

