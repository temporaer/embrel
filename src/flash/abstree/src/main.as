package src 
{
	/*
	 * Sprites sind so grafische Objekte, die andere enthalten können und solches Zeuch
	 */
	import flash.display.Sprite;
	import flash.events.Event;
	import src.Objs.xPoint;
	/*
	 * XML-Loader-Objekt
	 */
	import src.Loaders.xmlLoader;
	
	/**
	 * ...
	 * @author ...
	 */
	public class main extends Sprite
	{		
		private const pointURL:String = "data/points.xml";
		
		private var xml:xmlLoader = new xmlLoader();
		private var data:XML;
		
		
		public function main() 
		{
			this.xml.addEventListener( xmlLoader.XML_LOADED, onXMLLoaded, false, 0, true );
			this.xml.xmlurl = pointURL;			
		}
		
		private function onXMLLoaded(e:Event):void
		{
			this.data = this.xml.myXML;
			for ( var i:uint = 0; i < this.data.node.length(); i++ )
			{
				this.addChild( new xPoint( this.data.node[i] ) );
			}
		}
		
	}
	
}