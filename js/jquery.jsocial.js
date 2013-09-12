/*
Plugin Name: jSocial
Version: 1.0
Plugin URI: http://thatryan.com
Description: jQuery Plugin for Social Linking
Author: Ryan Olson
Author URI: http://thatryan.com

The beautiful icons used in this plugin are by Rogie King of http://www.komodomedia.com/
*/

/* PLUGIN USAGE */
/* 
	All that is needed to display the list is an empty div with the
	class of "social"
	
		<div class="social"></div>	

	At the bottom of your page, above the </body> tag, call the script like this, 
	
	<script type="text/javascript" charset="utf-8">
	$('.social').jsocial({

			//options go in here...
	});
	</script>
	
	The available options, (for now) are as follows, for most enter only your username
	unless otherwise specified. Be sure to surround with '' . For URL's, do not enter the begining http:// 
	
	twitter		: null, 	//add your user name
	facebook	: null, 	//add your profile URL
	flickr		: null, 	//add your user name
	friendfeed	: null, 	//add your user name
	delicious	: null, 	//add your user name
	digg		: null, 	//add your user name
	lastfm		: null, 	//add your user name
	linked		: null, 	//add your profile URL
	youtube		: null, 	//add your profile URL
	feed		: null,		//add your user name
	newPage		: null,		// true or false, opens link in new page
	center		: false,	//true or false, auto center the list
	inline		: false,	//true or false, display list horizontally
	small		: false		//true or false, use small (16x16)pixel icons
	
*/


(function($) {
	
	$.fn.jsocial = function(options) {

		var options = $.extend({
			
			twitter		: null, 
			facebook	: null, 
			flickr		: null, 
			friendfeed	: null, 
			delicious	: null, 
			digg		: null, 
			lastfm		: null, 
			linked		: null, 
			youtube		: null, 
			feed		: null,
			newPage		: null,
			center		: false,
			inline		: false,
			small		: false
			
		},options);
		
	
		return this.each(function() {
			
			if ( options.newPage ) {
				var target = 'target="_blank"';		//if user wants link to open in new page
			};
		
			var $this = $(this),		//cache 'this' , define variables
				total_links = 0,		//used to auto set the width of ul, for centering
				
				list = '<ul class = "social_links">';	//open the list element
						
if( options.twitter ) {list += '<li><a href="http://twitter.com/' + options.twitter + '" class="twitter" ' + target + '>Twitter</a></li>'; total_links ++;};

if( options.facebook ) {list += '<li><a href="http://' + options.facebook + '" class="facebook" ' + target + ' >Facebook</a></li>'; total_links ++;};

if( options.flickr ) {list += '<li><a href="http://flickr.com/' + options.flickr + '" class="flickr" ' + target + ' >Flickr</a></li>'; total_links ++;};

if( options.friendfeed ) {list += '<li><a href="http://friendfeed.com/' + options.friendfeed + '" class="friendfeed" ' + target + ' >FriendFeed</a></li>'; total_links ++;};

if( options.delicious ) {list += '<li><a href="http://delicious.com/' + options.delicious + '" class="delicious" ' + target + ' >Delicious</a></li>'; total_links ++;};

if( options.digg ) {list += '<li><a href="http://digg.com/users/' + options.digg + '" class="digg" ' + target + ' >Digg</a></li>'; total_links ++;};

if( options.lastfm ) {list += '<li><a href="http://www.last.fm/user/' + options.lastfm + '" class="lastfm" ' + target + ' >Last.FM</a></li>'; total_links ++;};

if( options.linked ) {list += '<li><a href="http://' + options.linked + '" class="linked-in" ' + target + ' >Linked-In</a></li>'; total_links ++;};

if( options.youtube ) {list += '<li><a href="http://' + options.youtube + '" class="youtube" ' + target + ' >YouTube</a></li>'; total_links ++;};

if( options.feed ) {list += '<li><a href="http://feeds.feedburner.com/' + options.feed + '" class="feed" ' + target + ' >Feedburner</a></li>'; total_links ++;};
			
				list += '</ul>';	//close the ul
			
			$this.append(list);	//add the list to the DOM
			
			if ( options.inline ) {
				$('ul.social_links li').css('float','left');	//add the float left to each list element
			};	//end inline if
			
			if ( options.small ) {
				$(this).removeClass('social').addClass('social_small');
			};	//end small icons if
			
			if ( options.center ) {		//if user wants social list centered
				
				if ( options.small) {	//if they chose small icons
					$('ul.social_links').css('width', total_links*26); //calculate the width of the ul for centering
					
				}
				else {		//else they chose large icons
				
				$('ul.social_links').css('width', total_links*42); //calculate the width of the ul for centering
				};
				
			};	//end centering if
			
		});  //end each
	
	} //end function
	
})(jQuery);


