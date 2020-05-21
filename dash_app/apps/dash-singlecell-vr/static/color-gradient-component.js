/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('color-gradient', {
    schema: {
        colors: {type: 'array'},
        height: {type: 'number'},
        width: {type: 'number'},
        verticalOffset: {type: 'number', default: 0},
    },

    init: function () {
        this.height = this.data.height;
        this.width = this.data.width;
        this.colors = this.data.colors;
        this.verticalOffset = -this.data.verticalOffset;
        const blocks = [];
        const blockHeight = 5 / this.colors.length;
        for (let i = 0; i < this.colors.length; i++) {
          const block = Utils.htmlToElement(`<a-entity position="0 ${this.verticalOffset + (blockHeight * i)} 0" geometry="primitive: plane; height: ${blockHeight}; width: ${this.width};" material="shader: flat; side:front; transparent: true; alphaTest: 0.5; color: ${this.colors[i]}"></a-entity>`);
          blocks.push(block);
        }
        this.el.setAttribute('geometry', `primitive: plane; height: ${this.height}; width: ${this.width};`);
        this.el.setAttribute('material', `shader: flat; side:front; transparent: true; opacity: 0; alphaTest: 0.5; color: white;`);  
        this.el.append(...blocks);
    }
});