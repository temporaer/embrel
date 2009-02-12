package src.Objs 
{
	import flash.display.Sprite;
	import flash.events.Event;
	import flash.events.MouseEvent;
	import flash.ui.Mouse;
	import flash.geom.ColorTransform;
	
	/**
	 * Punkt
	 * @author ...
	 */
	public class xPoint extends Sprite
	{
		private var node:XML;
		
		private var references:Array = new Array();
		private var color:Number = 0;
		private var id:uint = 0;
		private var size:uint = 0;
		
		
		public function xPoint(nodepoint:XML) 
		{
			this.node = nodepoint;
			this.addEventListener( Event.ADDED_TO_STAGE, draw, false, 0 , true );
			this.addEventListener( MouseEvent.MOUSE_OVER, onMouseOver, false, 0, true );
			this.addEventListener( MouseEvent.MOUSE_OUT, onMouseOut, false, 0, true );
		}
		
		
		private function draw( e:Event ):void
		{
			this.removeEventListener( Event.ADDED_TO_STAGE, draw );
			this.color = node.@color;
			this.name = node.@id;
			this.size = node.@size;
			this.x = node.@x;
			this.y = node.@y;
			for ( var i:uint = 0; i < node.child.length(); i++ )
			{
				this.references.push( node.child[i].@id );
			}
			this.graphics.beginFill( this.color, 1 );
			this.graphics.drawCircle( 0, 0, this.size * 0.5 );			
		}
		
		public function tint(color:Number = -1):void
		{		
			if ( color == -1 )
			{
				color = this.color;
			}
			
			var colorTransform:ColorTransform = this.transform.colorTransform;
			colorTransform.color = color;
			this.transform.colorTransform = colorTransform;		
			
		}
		
		
		private function onMouseOver(e:MouseEvent):void
		{
			for ( var i:uint = 0; i < this.references.length; i++)
			{				
				(this.parent.getChildByName( String(this.references[ i ] )) as xPoint).tint( this.color );
			}
		}
		
	
		
		private function onMouseOut(e:MouseEvent):void
		{
			for ( var i:uint = 0; i < this.references.length; i++)
			{
				(this.parent.getChildByName( this.references[ i ] ) as xPoint).tint();
			}
		}

		
		
		

		

		
	}
	
}