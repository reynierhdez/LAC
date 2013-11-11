//============================================================================================================================================================
// ==> BUSCA BIBLIOTECA E GOOGLE NO TOPO

/*function verificaBusca(frm){

	if(frm.tipo.value == "S")
	{
		frm.method = 'get';		
		frm.action = 'http://www.google.com/u/INPE';
	}
	else{
		frm.method = 'post';
		frm.action = '/resultado.php';
	}
	
}*/


//============================================================================================================================================================
// ==> BANNERS CBERS

function Banners()  //Funciona com refresh na página
{
    var MNews = new Array(); 
    MNews[0]= '<div style="background-color:#ffffff;"><a href="http://www.dgi.inpe.br:80/CDSR/" rel="externo"><img src="/imagens/cbers/brasilia.jpg" alt="" /></a><div style="background-color:#ffffff; padding:5px 5px 10px 5px; font-size:10px; font-weight:bold;">Bras&iacute;lia (DF)<br />Imagens CBERS-2 - C&acirc;mera CCD</div></div>';
    MNews[1]= '<div style="background-color:#ffffff;"><a href="http://www.dgi.inpe.br:80/CDSR/" rel="externo"><img src="/imagens/cbers/aracaju.jpg" alt="" /></a><div style="background-color:#ffffff; padding:5px 5px 10px 5px; font-size:10px; font-weight:bold;">Aracaju (SE)<br />Imagens CBERS-2 - C&acirc;mera CCD</div></div>';
    MNews[2]= '<div style="background-color:#ffffff;"><a href="http://www.dgi.inpe.br:80/CDSR/" rel="externo"><img src="/imagens/cbers/belem.jpg" alt="" /></a><div style="background-color:#ffffff; padding:5px 5px 10px 5px; font-size:10px; font-weight:bold;">Bel&eacute;m (PA)<br />Imagens CBERS-2 - C&acirc;mera CCD</div></div>';
	MNews[3]= '<div style="background-color:#ffffff;"><a href="http://www.dgi.inpe.br:80/CDSR/" rel="externo"><img src="/imagens/cbers/belo_horizonte.jpg" alt="" /></a><div style="background-color:#ffffff; padding:5px 5px 10px 5px; font-size:10px; font-weight:bold;">Belo Horizonte (MG)<br />Imagens CBERS-2 - C&acirc;mera CCD</div></div>';
	MNews[4]= '<div style="background-color:#ffffff;"><a href="http://www.dgi.inpe.br:80/CDSR/" rel="externo"><img src="/imagens/cbers/boa_vista.jpg" alt="" /></a><div style="background-color:#ffffff; padding:5px 5px 10px 5px; font-size:10px; font-weight:bold;">Boa Vista (RR)<br />Imagens CBERS-2 - C&acirc;mera CCD</div></div>';

      
   var Numero = Math.floor(Math.random()*5);
  
   document.write(MNews[Numero]);
}


//============================================================================================================================================================
// ==> AUMENTAR E DIMINUIR LETRA DO TEXTO
var tagAlvo = new Array('table');
var tamanhos = new Array( '8px','10px','11px','12px','14px','18px','22px' );
var tamanhoInicial = 2;
 
function mudaTamanho( idAlvo, acao ){
  if (!document.getElementById) return
  var selecionados = null,tamanho = tamanhoInicial,i,j,tagsAlvo;
  tamanho += acao;
  if ( acao == 0 ) tamanho = 2;
  if ( tamanho < 0 ) tamanho = 0;
  if ( tamanho > 6 ) tamanho = 6;
  tamanhoInicial = tamanho;
  if ( !( selecionados = document.getElementById( idAlvo ) ) ) selecionados = document.getElementsByTagName( idAlvo )[ 0 ];
  
  selecionados.style.fontSize = tamanhos[ tamanho ];
  
  for ( i = 0; i < tagAlvo.length; i++ ){
    tagsAlvo = selecionados.getElementsByTagName( tagAlvo[ i ] );
    for ( j = 0; j < tagsAlvo.length; j++ ) tagsAlvo[ j ].style.fontSize = tamanhos[ tamanho ];
  }
}

//==========================================================================================================================================================
// ==> ABRIR SUBMENU
startList = function() 
{
	/*// ==> CRIAR SUBMENU
	var navItems = document.getElementById("menu_principal").getElementsByTagName('li');
	for (var i=0; i< navItems.length; i++) {
		if(navItems[i].className == "submenu") {
			var subnavItem = navItems[i].getElementsByTagName('ul')[0];
			if(navItems[i].getElementsByTagName('ul')[0] != null) {
				navItems[i].onmouseover = function() {
					this.getElementsByTagName('ul')[0].style.display = "block";
				}
				navItems[i].onmouseout = function() {
					this.getElementsByTagName('ul')[0].style.display = "none";
				}
			}
		}
	}*/
	
	// CRIAR LINKS EXTERNOS
	if(document.getElementsByTagName) {
		var anchors = document.getElementsByTagName('a');
		
		for(var i=0; i<anchors.length; i++) {
			var anchor = anchors[i];
			if(anchor.getAttribute("href") && anchor.getAttribute('rel')=="externo") { // <-- é necessário inserir rel="externo" no link
				anchor.target = '_blank';
				var title = anchor.title + ''; // <-- Insere este texto no final do Title do link
				anchor.title = title;
			}		
			
		}
		
		var anchors2 = document.getElementsByTagName('form');
		
		for(var i=0; i<anchors2.length; i++) {
			var anchor = anchors2[i];
			if(anchor.getAttribute('class')=="externo") { // <-- é necessário inserir class="externo" no link
				anchor.target = '_blank';
				var title = anchor.title + ''; // <-- Insere este texto no final do Title do link
				anchor.title = title;
			}		
			
		}

		var inputs = document.getElementsByTagName('input');
		if (inputs.length > 0) {
			inputs[0].focus();
		}
	}

}
window.onload=startList; 

//========================================================================================================================================================
// ==> MUDAR IMAGEM DE COR. (LEIA MAIS)
function mudar(qual, i)
{
	var img = document.getElementById("mais"+i);
	img.src = qual;
}

//=========================================================================================================================================================
// ==> CATÁLAGO DE IMAGENS (CBERS - LANDSAT - GOES - MSG)
function imagens(qual, id, link)
{
	var img = document.getElementById("imagem");
	img.src = qual;
	
	document.getElementById("n1").style.background = "#e6e3e5";
	document.getElementById("n2").style.background = "#e6e3e5";
	document.getElementById("n3").style.background = "#e6e3e5";
	document.getElementById("n4").style.background = "#e6e3e5";
	
	document.getElementById(id).style.background = "#c4c2c3";
	
	document.getElementById("link").href=link;
}

//======================================================================================================================================================
// ==> COLOCA RELÓGIO NO TOPO DA PÁGINA
/*var req;

function loadXMLDoc(url)
{
    req = null;
    // Procura por um objeto nativo (Mozilla/Safari)
    if (window.XMLHttpRequest) {
        req = new XMLHttpRequest();
        req.onreadystatechange = processReqChange;
        req.open("GET", url, true);
        req.send(null);
    // Procura por uma versão ActiveX (IE)
    } else if (window.ActiveXObject) {
        try {
req = new ActiveXObject("Msxml2.XMLHTTP.4.0"); //alert(req);
} catch(e) {
try {
req = new ActiveXObject("Msxml2.XMLHTTP.3.0"); //alert(req);
} catch(e) {
try {
req = new ActiveXObject("Msxml2.XMLHTTP"); //alert(req);
} catch(e) {
try {
req = new ActiveXObject("Microsoft.XMLHTTP"); //alert(req);
} catch(e) {
req = false;
}
}
}
}
        if (req) {
            req.onreadystatechange = processReqChange;
            req.open("GET", url, true);
            req.send();
        }
    }
}

function processReqChange()
{
    // apenas quando o estado for "completado"
    if (req.readyState == 4) {
        // apenas se o servidor retornar "OK"
        if (req.status == 200) {
            // procura pela div id="news" e insere o conteudo
            // retornado nela, como texto HTML
            document.getElementById('relogio').innerHTML = req.responseText;
        } 
		//else {
        //    alert("Erro na hora:n" + req.statusText);
        //}
    }
}

function buscarTempo()
{
    loadXMLDoc("/include/relogio.php");
}

// Recarrega a cada 60000 milissegundo (60 segundos)
setInterval("buscarTempo()", 1000); 
*/

//======================================================================================================================================================
// ==> ANÁLISE PELO GOOGLE ANALYTICS
//var _gaq = _gaq || [];
//_gaq.push(['_setAccount', 'UA-5057532-1']);
//_gaq.push(['_trackPageview']);

//  (function() {
//    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
//    ga.src = ('https:' == document.location.protocol ? 'https://ssl' :
//'http://www') + '.google-analytics.com/ga.js';
//    var s = document.getElementsByTagName('script')[0];
//s.parentNode.insertBefore(ga, s);
//  })();


