// JavaScript Document
$(document).on("ready",configurarApp);

function configurarApp()
{
var canvas = document.getElementById('myCanvas');
var context = canvas.getContext('2d');
var x = 80;
var y = 110;
context.font = '41pt TwCen';
context.lineWidth = 1;
// stroke color
context.fillStyle = 'white';
context.fillText('Computação Aplicada', x, y);
context.strokeStyle = 'black';
context.strokeText('Computação Aplicada', x, y);
}